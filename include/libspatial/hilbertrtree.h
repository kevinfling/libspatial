/**
 * @file   hilbertrtree.h
 * @brief  Hilbert R-Tree with R* splitting, forced reinsert, and bulk loading
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * Hilbert R-Tree combining:
 * - Hilbert curve ordering for cache-friendly spatial locality
 * - R* split heuristic (minimum margin/overlap/area)
 * - Forced reinsertion on overflow (R* algorithm)
 * - SoA node layout for vectorization-friendly traversal
 * - Radix sort for O(n) Hilbert ordering on bulk load
 * - Stack-based (non-recursive) search
 * - Standard spatial_allocator / arena / pool / memento integration
 */

#ifndef SPATIAL_HILBERTRTREE_H
#define SPATIAL_HILBERTRTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

/* ============================================================================
 * Configuration
 * ============================================================================ */

#ifndef SPATIAL_HILBERTRTREE_DIMS
    #define SPATIAL_HILBERTRTREE_DIMS 2
#endif

#ifndef SPATIAL_HILBERTRTREE_MAX_CHILDREN
    #define SPATIAL_HILBERTRTREE_MAX_CHILDREN 16
#endif

#ifndef SPATIAL_HILBERTRTREE_MIN_CHILDREN
    #define SPATIAL_HILBERTRTREE_MIN_CHILDREN (SPATIAL_HILBERTRTREE_MAX_CHILDREN / 4)
#endif

#ifndef SPATIAL_HILBERTRTREE_REINSERT_PCT
    #define SPATIAL_HILBERTRTREE_REINSERT_PCT 30
#endif

#ifndef SPATIAL_HILBERTRTREE_MAX_REINSERT_DEPTH
    #define SPATIAL_HILBERTRTREE_MAX_REINSERT_DEPTH 8
#endif

SPATIAL_STATIC_ASSERT(hrt_dims_valid,
    SPATIAL_HILBERTRTREE_DIMS >= 2 && SPATIAL_HILBERTRTREE_DIMS <= 16);
SPATIAL_STATIC_ASSERT(hrt_max_children_valid,
    SPATIAL_HILBERTRTREE_MAX_CHILDREN >= 4 && SPATIAL_HILBERTRTREE_MAX_CHILDREN <= 256);

/* ============================================================================
 * Internal types
 * ============================================================================ */

/* Forward declaration */
struct spatial_hilbertrtree_node;

/* Temporary entry used during splits and reinserts */
typedef struct {
    spatial_num_t min[SPATIAL_HILBERTRTREE_DIMS];
    spatial_num_t max[SPATIAL_HILBERTRTREE_DIMS];
    uint64_t hilbert;
    union {
        struct spatial_hilbertrtree_node *child;
        spatial_data_t data;
    } ptr;
    bool is_leaf; /* true = data entry, false = child pointer */
} spatial_hrt_entry;

/*
 * Node uses Structure-of-Arrays layout so the hot traversal arrays
 * (mins, maxs, hilberts) sit in contiguous memory for vectorization.
 * children/datas are accessed only when a box test passes.
 */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_hilbertrtree_node {
    spatial_num_t mins[SPATIAL_HILBERTRTREE_MAX_CHILDREN][SPATIAL_HILBERTRTREE_DIMS];
    spatial_num_t maxs[SPATIAL_HILBERTRTREE_MAX_CHILDREN][SPATIAL_HILBERTRTREE_DIMS];
    uint64_t hilberts[SPATIAL_HILBERTRTREE_MAX_CHILDREN];

    struct spatial_hilbertrtree_node *children[SPATIAL_HILBERTRTREE_MAX_CHILDREN];
    spatial_data_t datas[SPATIAL_HILBERTRTREE_MAX_CHILDREN];

    /* Aggregated bounds for quick node-level rejection */
    spatial_num_t node_min[SPATIAL_HILBERTRTREE_DIMS];
    spatial_num_t node_max[SPATIAL_HILBERTRTREE_DIMS];

    uint32_t count;
    uint32_t level;
    bool is_leaf;
} spatial_hilbertrtree_node;

/* Main tree structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_hilbertrtree {
    spatial_hilbertrtree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;

    /* Global bounds for Hilbert normalization */
    spatial_num_t global_min[SPATIAL_HILBERTRTREE_DIMS];
    spatial_num_t global_max[SPATIAL_HILBERTRTREE_DIMS];
    bool bounds_valid;

    int count;
    int height;
    bool relaxed_atomics;
} spatial_hilbertrtree;

/* Iterator callback — return true to continue, false to stop early */
typedef bool (*spatial_hilbertrtree_iter_fn)(const spatial_num_t *min,
                                              const spatial_num_t *max,
                                              spatial_data_t data,
                                              void *udata);

/* Forward declaration for public API */
SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_new_with_allocator(
    const spatial_allocator *alloc);

/* ============================================================================
 * Hilbert curve encoding
 * ============================================================================ */

/*
 * Encodes up to SPATIAL_HILBERTRTREE_DIMS coordinates as a 64-bit Hilbert
 * index by normalising to [0,1], applying Gray code, and interleaving bits.
 * BITS_PER_DIM = 63/dims so the result always fits in uint64_t.
 */
SPATIAL_INLINE uint64_t spatial_hrt_hilbert_encode(
    const spatial_num_t *SPATIAL_RESTRICT coords,
    const spatial_num_t *SPATIAL_RESTRICT gmin,
    const spatial_num_t *SPATIAL_RESTRICT gmax,
    int dims)
{
    const int BITS = 63 / dims;
    const uint64_t MAXVAL = (1ULL << BITS) - 1ULL;
    uint64_t h = 0;

    for (int d = 0; d < dims; d++) {
        spatial_num_t range = gmax[d] - gmin[d];
        spatial_num_t norm = (range > (spatial_num_t)0.0) ?
            (coords[d] - gmin[d]) / range : (spatial_num_t)0.5;

        if (norm < (spatial_num_t)0.0) norm = (spatial_num_t)0.0;
        if (norm > (spatial_num_t)1.0) norm = (spatial_num_t)1.0;

        uint64_t ci = (uint64_t)(norm * (spatial_num_t)MAXVAL);
        uint64_t gray = ci ^ (ci >> 1);

        for (int b = 0; b < BITS; b++) {
            h |= ((gray >> b) & 1ULL) << (d + b * dims);
        }
    }

    return h;
}

SPATIAL_INLINE uint64_t spatial_hrt_hilbert_for_tree(
    const spatial_hilbertrtree *SPATIAL_RESTRICT rt,
    const spatial_num_t *SPATIAL_RESTRICT min,
    const spatial_num_t *SPATIAL_RESTRICT max)
{
    spatial_num_t center[SPATIAL_HILBERTRTREE_DIMS];
    for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
        center[d] = (min[d] + max[d]) * (spatial_num_t)0.5;
    }

    if (rt->bounds_valid) {
        return spatial_hrt_hilbert_encode(center, rt->global_min, rt->global_max,
                                          SPATIAL_HILBERTRTREE_DIMS);
    }

    /* Unit cube fallback when bounds not yet set */
    spatial_num_t unit_min[SPATIAL_HILBERTRTREE_DIMS];
    spatial_num_t unit_max[SPATIAL_HILBERTRTREE_DIMS];
    for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
        unit_min[d] = (spatial_num_t)0.0;
        unit_max[d] = (spatial_num_t)1.0;
    }
    return spatial_hrt_hilbert_encode(center, unit_min, unit_max,
                                      SPATIAL_HILBERTRTREE_DIMS);
}

/* ============================================================================
 * Node management
 * ============================================================================ */

SPATIAL_INLINE void spatial_hrt_node_update_bounds(
    spatial_hilbertrtree_node *SPATIAL_RESTRICT node)
{
    if (node->count == 0) {
        spatial_bbox_init(node->node_min, node->node_max, SPATIAL_HILBERTRTREE_DIMS);
        return;
    }
    for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
        node->node_min[d] = node->mins[0][d];
        node->node_max[d] = node->maxs[0][d];
    }
    for (uint32_t i = 1; i < node->count; i++) {
        for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
            node->node_min[d] = spatial_min(node->node_min[d], node->mins[i][d]);
            node->node_max[d] = spatial_max(node->node_max[d], node->maxs[i][d]);
        }
    }
}

SPATIAL_INLINE spatial_hilbertrtree_node *spatial_hrt_node_new(
    bool is_leaf, uint32_t level, const spatial_allocator *SPATIAL_RESTRICT alloc)
{
    spatial_hilbertrtree_node *node = (spatial_hilbertrtree_node *)alloc->malloc(
        sizeof(spatial_hilbertrtree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;

    memset(node, 0, sizeof(spatial_hilbertrtree_node));
    node->level   = level;
    node->is_leaf = is_leaf;
    spatial_bbox_init(node->node_min, node->node_max, SPATIAL_HILBERTRTREE_DIMS);

    return node;
}

SPATIAL_INLINE void spatial_hrt_node_free(
    spatial_hilbertrtree_node *node,
    const spatial_allocator *alloc,
    const spatial_item_callbacks *cb,
    void *udata);

SPATIAL_INLINE void spatial_hrt_node_free(
    spatial_hilbertrtree_node *node,
    const spatial_allocator *alloc,
    const spatial_item_callbacks *cb,
    void *udata)
{
    if (SPATIAL_UNLIKELY(!node)) return;

    if (node->is_leaf) {
        if (cb && cb->free) {
            for (uint32_t i = 0; i < node->count; i++) {
                cb->free(node->datas[i], udata);
            }
        }
    } else {
        for (uint32_t i = 0; i < node->count; i++) {
            spatial_hrt_node_free(node->children[i], alloc, cb, udata);
        }
    }

    alloc->free(node, alloc->udata);
}

/* ============================================================================
 * R* split algorithm
 * ============================================================================ */

/*
 * Compute margin (sum of side lengths) of a bounding box — used by
 * ChooseSplitAxis to pick the axis that minimises perimeter.
 */
SPATIAL_INLINE spatial_num_t spatial_hrt_bbox_margin(
    const spatial_num_t *SPATIAL_RESTRICT mn,
    const spatial_num_t *SPATIAL_RESTRICT mx)
{
    spatial_num_t m = (spatial_num_t)0.0;
    for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
        m += mx[d] - mn[d];
    }
    return (spatial_num_t)2.0 * m;
}

SPATIAL_INLINE spatial_num_t spatial_hrt_bbox_overlap(
    const spatial_num_t *SPATIAL_RESTRICT a_min, const spatial_num_t *SPATIAL_RESTRICT a_max,
    const spatial_num_t *SPATIAL_RESTRICT b_min, const spatial_num_t *SPATIAL_RESTRICT b_max)
{
    spatial_num_t ov = (spatial_num_t)1.0;
    for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
        spatial_num_t o = spatial_min(a_max[d], b_max[d]) - spatial_max(a_min[d], b_min[d]);
        if (o <= (spatial_num_t)0.0) return (spatial_num_t)0.0;
        ov *= o;
    }
    return ov;
}

/*
 * Sort entries[] by their min[axis] using insertion sort.
 * HRT_MAX_CHILDREN <= 256, so insertion sort is cache-friendly here.
 */
SPATIAL_INLINE void spatial_hrt_sort_by_axis(
    spatial_hrt_entry *SPATIAL_RESTRICT entries, int count, int axis)
{
    for (int i = 1; i < count; i++) {
        spatial_hrt_entry key = entries[i];
        int j = i - 1;
        while (j >= 0 && entries[j].min[axis] > key.min[axis]) {
            entries[j + 1] = entries[j];
            j--;
        }
        entries[j + 1] = key;
    }
}

/*
 * R* ChooseSplitAxis: returns the axis whose sorted distributions have
 * the minimum total margin sum.
 */
SPATIAL_INLINE int spatial_hrt_choose_split_axis(
    spatial_hrt_entry *SPATIAL_RESTRICT entries, int count)
{
    spatial_num_t best_margin = (spatial_num_t)DBL_MAX;
    int best_axis = 0;

    for (int axis = 0; axis < SPATIAL_HILBERTRTREE_DIMS; axis++) {
        spatial_hrt_sort_by_axis(entries, count, axis);

        spatial_num_t margin_sum = (spatial_num_t)0.0;
        for (int k = SPATIAL_HILBERTRTREE_MIN_CHILDREN;
             k <= count - SPATIAL_HILBERTRTREE_MIN_CHILDREN; k++) {

            /* Group 1 bounds [0, k) */
            spatial_num_t g1_min[SPATIAL_HILBERTRTREE_DIMS];
            spatial_num_t g1_max[SPATIAL_HILBERTRTREE_DIMS];
            memcpy(g1_min, entries[0].min, sizeof(g1_min));
            memcpy(g1_max, entries[0].max, sizeof(g1_max));
            for (int i = 1; i < k; i++) {
                for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                    g1_min[d] = spatial_min(g1_min[d], entries[i].min[d]);
                    g1_max[d] = spatial_max(g1_max[d], entries[i].max[d]);
                }
            }

            /* Group 2 bounds [k, count) */
            spatial_num_t g2_min[SPATIAL_HILBERTRTREE_DIMS];
            spatial_num_t g2_max[SPATIAL_HILBERTRTREE_DIMS];
            memcpy(g2_min, entries[k].min, sizeof(g2_min));
            memcpy(g2_max, entries[k].max, sizeof(g2_max));
            for (int i = k + 1; i < count; i++) {
                for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                    g2_min[d] = spatial_min(g2_min[d], entries[i].min[d]);
                    g2_max[d] = spatial_max(g2_max[d], entries[i].max[d]);
                }
            }

            margin_sum += spatial_hrt_bbox_margin(g1_min, g1_max)
                        + spatial_hrt_bbox_margin(g2_min, g2_max);
        }

        if (margin_sum < best_margin) {
            best_margin = margin_sum;
            best_axis   = axis;
        }
    }

    return best_axis;
}

/*
 * R* ChooseSplitIndex: given entries sorted on best_axis, find the split
 * point that minimises overlap, breaking ties by minimum total area.
 */
SPATIAL_INLINE int spatial_hrt_choose_split_index(
    spatial_hrt_entry *SPATIAL_RESTRICT entries, int count, int axis)
{
    spatial_hrt_sort_by_axis(entries, count, axis);

    spatial_num_t best_overlap = (spatial_num_t)DBL_MAX;
    spatial_num_t best_area    = (spatial_num_t)DBL_MAX;
    int best_k = SPATIAL_HILBERTRTREE_MIN_CHILDREN;

    for (int k = SPATIAL_HILBERTRTREE_MIN_CHILDREN;
         k <= count - SPATIAL_HILBERTRTREE_MIN_CHILDREN; k++) {

        spatial_num_t g1_min[SPATIAL_HILBERTRTREE_DIMS];
        spatial_num_t g1_max[SPATIAL_HILBERTRTREE_DIMS];
        memcpy(g1_min, entries[0].min, sizeof(g1_min));
        memcpy(g1_max, entries[0].max, sizeof(g1_max));
        for (int i = 1; i < k; i++) {
            for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                g1_min[d] = spatial_min(g1_min[d], entries[i].min[d]);
                g1_max[d] = spatial_max(g1_max[d], entries[i].max[d]);
            }
        }

        spatial_num_t g2_min[SPATIAL_HILBERTRTREE_DIMS];
        spatial_num_t g2_max[SPATIAL_HILBERTRTREE_DIMS];
        memcpy(g2_min, entries[k].min, sizeof(g2_min));
        memcpy(g2_max, entries[k].max, sizeof(g2_max));
        for (int i = k + 1; i < count; i++) {
            for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                g2_min[d] = spatial_min(g2_min[d], entries[i].min[d]);
                g2_max[d] = spatial_max(g2_max[d], entries[i].max[d]);
            }
        }

        spatial_num_t overlap = spatial_hrt_bbox_overlap(g1_min, g1_max, g2_min, g2_max);
        spatial_num_t area    = spatial_bbox_area(g1_min, g1_max, SPATIAL_HILBERTRTREE_DIMS)
                              + spatial_bbox_area(g2_min, g2_max, SPATIAL_HILBERTRTREE_DIMS);

        if (overlap < best_overlap ||
            (!(overlap > best_overlap) && area < best_area)) {
            best_overlap = overlap;
            best_area    = area;
            best_k       = k;
        }
    }

    return best_k;
}

/*
 * Unpack a node's SoA arrays into a flat entry array for split processing,
 * run R* split, then repack into the original node and a new sibling.
 */
SPATIAL_INLINE spatial_hilbertrtree_node *spatial_hrt_split_node(
    spatial_hilbertrtree_node *SPATIAL_RESTRICT node,
    const spatial_allocator *SPATIAL_RESTRICT alloc)
{
    int count = (int)node->count;

    /* Unpack SoA -> AoS */
    spatial_hrt_entry entries[SPATIAL_HILBERTRTREE_MAX_CHILDREN];
    for (int i = 0; i < count; i++) {
        memcpy(entries[i].min,  node->mins[i],   sizeof(entries[i].min));
        memcpy(entries[i].max,  node->maxs[i],   sizeof(entries[i].max));
        entries[i].hilbert   = node->hilberts[i];
        entries[i].is_leaf   = node->is_leaf;
        if (node->is_leaf) {
            entries[i].ptr.data  = node->datas[i];
        } else {
            entries[i].ptr.child = node->children[i];
        }
    }

    /* R* choose split */
    int axis      = spatial_hrt_choose_split_axis(entries, count);
    int split_idx = spatial_hrt_choose_split_index(entries, count, axis);

    /* Allocate sibling */
    spatial_hilbertrtree_node *sibling = spatial_hrt_node_new(
        node->is_leaf, node->level, alloc);
    if (SPATIAL_UNLIKELY(!sibling)) return NULL;

    /* Repack first group back into original node */
    node->count = (uint32_t)split_idx;
    for (int i = 0; i < split_idx; i++) {
        memcpy(node->mins[i],  entries[i].min,  sizeof(node->mins[i]));
        memcpy(node->maxs[i],  entries[i].max,  sizeof(node->maxs[i]));
        node->hilberts[i] = entries[i].hilbert;
        if (node->is_leaf) {
            node->datas[i]    = entries[i].ptr.data;
        } else {
            node->children[i] = entries[i].ptr.child;
        }
    }
    spatial_hrt_node_update_bounds(node);

    /* Pack second group into sibling */
    int sib_count = count - split_idx;
    sibling->count = (uint32_t)sib_count;
    for (int i = 0; i < sib_count; i++) {
        memcpy(sibling->mins[i],  entries[split_idx + i].min,  sizeof(sibling->mins[i]));
        memcpy(sibling->maxs[i],  entries[split_idx + i].max,  sizeof(sibling->maxs[i]));
        sibling->hilberts[i] = entries[split_idx + i].hilbert;
        if (sibling->is_leaf) {
            sibling->datas[i]    = entries[split_idx + i].ptr.data;
        } else {
            sibling->children[i] = entries[split_idx + i].ptr.child;
        }
    }
    spatial_hrt_node_update_bounds(sibling);

    return sibling;
}

/* ============================================================================
 * Path-recording descent + iterative upward insert / split propagation
 * ============================================================================ */

/*
 * Path recorded during descent so we can propagate splits upward without
 * re-descending from root.  Max height = 64 levels is more than enough for
 * any tree with MAX_CHILDREN = 16 and billions of items.
 */
#define SPATIAL_HRT_MAX_PATH 64

typedef struct {
    spatial_hrt_entry entry;
    int level;
} spatial_hrt_reinsert_item;

typedef struct {
    spatial_hrt_reinsert_item items[SPATIAL_HILBERTRTREE_MAX_CHILDREN];
    int count;
    int depth;
} spatial_hrt_reinsert_queue;

/*
 * Insert `entry` at `target_level` (0 = leaf level) into the subtree rooted
 * at `path[path_len-1]`, propagating splits upward through the path.
 *
 * path[0] is root, path[path_len-1] is the target node.
 *
 * Returns true on success.  On split of the root a new root is installed in
 * `rt` and rt->height is incremented.
 */
static bool spatial_hrt_insert_at(
    spatial_hilbertrtree *rt,
    spatial_hilbertrtree_node **path,
    int path_len,
    const spatial_hrt_entry *entry,
    bool is_leaf_entry,
    spatial_hrt_reinsert_queue *queue)
{
    /* Work from the target node (bottom of path) upward */
    spatial_hrt_entry to_insert = *entry;
    to_insert.is_leaf = is_leaf_entry;

    for (int pi = path_len - 1; pi >= 0; pi--) {
        spatial_hilbertrtree_node *node = path[pi];

        /* ---------------------------------------------------------------
         * Fast path: node has room — insert and propagate bounds upward
         * --------------------------------------------------------------- */
        if (node->count < SPATIAL_HILBERTRTREE_MAX_CHILDREN) {
            uint32_t idx = node->count;
            memcpy(node->mins[idx],    to_insert.min, sizeof(node->mins[idx]));
            memcpy(node->maxs[idx],    to_insert.max, sizeof(node->maxs[idx]));
            node->hilberts[idx] = to_insert.hilbert;
            if (to_insert.is_leaf) {
                node->datas[idx]    = to_insert.ptr.data;
            } else {
                node->children[idx] = to_insert.ptr.child;
            }
            node->count++;

            /* Incremental bounds update for this node */
            for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                node->node_min[d] = spatial_min(node->node_min[d], to_insert.min[d]);
                node->node_max[d] = spatial_max(node->node_max[d], to_insert.max[d]);
            }

            /* Propagate updated bounds to ancestors */
            for (int ai = pi - 1; ai >= 0; ai--) {
                spatial_hilbertrtree_node *anc = path[ai];
                /* Find which child slot holds path[ai+1] and update its bbox */
                spatial_hilbertrtree_node *child = path[ai + 1];
                for (uint32_t ci = 0; ci < anc->count; ci++) {
                    if (anc->children[ci] == child) {
                        memcpy(anc->mins[ci], child->node_min, sizeof(anc->mins[ci]));
                        memcpy(anc->maxs[ci], child->node_max, sizeof(anc->maxs[ci]));
                        anc->hilberts[ci] = spatial_hrt_hilbert_for_tree(
                            rt, child->node_min, child->node_max);
                        /* Recompute ancestor's own bounds */
                        spatial_hrt_node_update_bounds(anc);
                        break;
                    }
                }
            }
            return true;
        }

        /* ---------------------------------------------------------------
         * Node full at leaf level — R* forced reinsert (first overflow only)
         * --------------------------------------------------------------- */
        if (node->is_leaf && queue->depth < SPATIAL_HILBERTRTREE_MAX_REINSERT_DEPTH) {
            int reinsert_n = ((int)node->count * SPATIAL_HILBERTRTREE_REINSERT_PCT) / 100;
            if (reinsert_n < 1) reinsert_n = 1;
            if (reinsert_n > (int)node->count - SPATIAL_HILBERTRTREE_MIN_CHILDREN) {
                reinsert_n = (int)node->count - SPATIAL_HILBERTRTREE_MIN_CHILDREN;
            }

            /* Distance of each entry's center from node center */
            spatial_num_t nc[SPATIAL_HILBERTRTREE_DIMS];
            for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                nc[d] = (node->node_min[d] + node->node_max[d]) * (spatial_num_t)0.5;
            }

            typedef struct { int idx; spatial_num_t dist_sq; } dist_rec;
            dist_rec dists[SPATIAL_HILBERTRTREE_MAX_CHILDREN];
            for (uint32_t i = 0; i < node->count; i++) {
                spatial_num_t dsq = (spatial_num_t)0.0;
                for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                    spatial_num_t ec = (node->mins[i][d] + node->maxs[i][d]) * (spatial_num_t)0.5;
                    spatial_num_t diff = ec - nc[d];
                    dsq += diff * diff;
                }
                dists[i].idx     = (int)i;
                dists[i].dist_sq = dsq;
            }

            /* Insertion sort: largest dist_sq first */
            for (int i = 1; i < (int)node->count; i++) {
                dist_rec key = dists[i];
                int j = i - 1;
                while (j >= 0 && dists[j].dist_sq < key.dist_sq) {
                    dists[j + 1] = dists[j];
                    j--;
                }
                dists[j + 1] = key;
            }

            /* Queue furthest entries for deferred reinsert */
            bool removed[SPATIAL_HILBERTRTREE_MAX_CHILDREN];
            memset(removed, 0, sizeof(removed));
            for (int i = 0; i < reinsert_n && queue->count < SPATIAL_HILBERTRTREE_MAX_CHILDREN; i++) {
                int idx = dists[i].idx;
                spatial_hrt_reinsert_item *item = &queue->items[queue->count++];
                memcpy(item->entry.min,  node->mins[idx],  sizeof(item->entry.min));
                memcpy(item->entry.max,  node->maxs[idx],  sizeof(item->entry.max));
                item->entry.hilbert    = node->hilberts[idx];
                item->entry.ptr.data   = node->datas[idx];
                item->entry.is_leaf    = true;
                item->level            = 0;
                removed[idx] = true;
            }

            /* Compact node */
            uint32_t new_count = 0;
            for (uint32_t i = 0; i < node->count; i++) {
                if (!removed[i]) {
                    if (new_count != i) {
                        memcpy(node->mins[new_count],  node->mins[i],  sizeof(node->mins[0]));
                        memcpy(node->maxs[new_count],  node->maxs[i],  sizeof(node->maxs[0]));
                        node->hilberts[new_count] = node->hilberts[i];
                        node->datas[new_count]    = node->datas[i];
                    }
                    new_count++;
                }
            }
            node->count = new_count;
            spatial_hrt_node_update_bounds(node);

            /* Retry inserting the original entry into the now-uncrowded node */
            queue->depth++;
            bool ok = spatial_hrt_insert_at(rt, path, pi + 1, &to_insert,
                                             to_insert.is_leaf, queue);
            queue->depth--;
            return ok;
        }

        /* ---------------------------------------------------------------
         * Node full — split it
         * --------------------------------------------------------------- */
        spatial_hilbertrtree_node *sibling = spatial_hrt_split_node(node, &rt->alloc);
        if (SPATIAL_UNLIKELY(!sibling)) return false;

        /* Place the new entry into whichever half it belongs to by Hilbert order */
        spatial_hilbertrtree_node *target = node;
        if (node->count > 0 && sibling->count > 0) {
            uint64_t last_h = node->hilberts[node->count - 1];
            uint64_t sib_h  = sibling->hilberts[0];
            uint64_t d1 = (last_h > to_insert.hilbert) ?
                (last_h - to_insert.hilbert) : (to_insert.hilbert - last_h);
            uint64_t d2 = (sib_h > to_insert.hilbert) ?
                (sib_h - to_insert.hilbert) : (to_insert.hilbert - sib_h);
            if (d2 < d1) target = sibling;
        }

        uint32_t sidx = target->count;
        memcpy(target->mins[sidx],  to_insert.min, sizeof(target->mins[sidx]));
        memcpy(target->maxs[sidx],  to_insert.max, sizeof(target->maxs[sidx]));
        target->hilberts[sidx] = to_insert.hilbert;
        if (to_insert.is_leaf) {
            target->datas[sidx]    = to_insert.ptr.data;
        } else {
            target->children[sidx] = to_insert.ptr.child;
        }
        target->count++;
        spatial_hrt_node_update_bounds(target);
        spatial_hrt_node_update_bounds(node); /* may have lost entries */

        /* Are we at the root? */
        if (pi == 0) {
            /* Grow tree: create new root with node + sibling as children */
            spatial_hilbertrtree_node *new_root =
                spatial_hrt_node_new(false, (uint32_t)rt->height + 1, &rt->alloc);
            if (SPATIAL_UNLIKELY(!new_root)) {
                spatial_hrt_node_free(sibling, &rt->alloc, NULL, NULL);
                return false;
            }

            memcpy(new_root->mins[0],  node->node_min, sizeof(new_root->mins[0]));
            memcpy(new_root->maxs[0],  node->node_max, sizeof(new_root->maxs[0]));
            new_root->hilberts[0] = spatial_hrt_hilbert_for_tree(rt, node->node_min, node->node_max);
            new_root->children[0] = node;

            memcpy(new_root->mins[1],  sibling->node_min, sizeof(new_root->mins[1]));
            memcpy(new_root->maxs[1],  sibling->node_max, sizeof(new_root->maxs[1]));
            new_root->hilberts[1] = spatial_hrt_hilbert_for_tree(rt, sibling->node_min, sibling->node_max);
            new_root->children[1] = sibling;

            new_root->count = 2;
            spatial_hrt_node_update_bounds(new_root);

            rt->root   = new_root;
            rt->height++;
            return true;
        }

        /*
         * Not the root: prepare a sibling entry to insert into the parent
         * (path[pi-1]).  Loop back up to pi-1.
         */
        memcpy(to_insert.min,  sibling->node_min, sizeof(to_insert.min));
        memcpy(to_insert.max,  sibling->node_max, sizeof(to_insert.max));
        to_insert.hilbert    = spatial_hrt_hilbert_for_tree(rt, sibling->node_min, sibling->node_max);
        to_insert.ptr.child  = sibling;
        to_insert.is_leaf    = false;

        /* Before continuing upward, update the parent's slot for `node`
         * so its cached bbox reflects any entries that moved to the sibling. */
        spatial_hilbertrtree_node *parent = path[pi - 1];
        for (uint32_t ci = 0; ci < parent->count; ci++) {
            if (parent->children[ci] == node) {
                memcpy(parent->mins[ci], node->node_min, sizeof(parent->mins[ci]));
                memcpy(parent->maxs[ci], node->node_max, sizeof(parent->maxs[ci]));
                parent->hilberts[ci] = spatial_hrt_hilbert_for_tree(
                    rt, node->node_min, node->node_max);
                break;
            }
        }
        /* Continue loop: pi-- will try to insert to_insert (the sibling)
         * into path[pi-1] (the parent). */
    }

    /* Should never reach here — root case is handled inside the loop */
    return false;
}

/*
 * Descend from root to the leaf level, recording the path, then call
 * spatial_hrt_insert_at to do the actual insert + split propagation.
 */
static bool spatial_hrt_insert_entry(
    spatial_hilbertrtree *rt,
    const spatial_hrt_entry *entry,
    spatial_hrt_reinsert_queue *queue)
{
    /* Record the path from root to target leaf */
    spatial_hilbertrtree_node *path[SPATIAL_HRT_MAX_PATH];
    int path_len = 0;

    spatial_hilbertrtree_node *node = rt->root;
    while (true) {
        if (SPATIAL_UNLIKELY(path_len >= SPATIAL_HRT_MAX_PATH)) return false;
        path[path_len++] = node;

        if (node->is_leaf) break;

        /* Descend to child with minimum Hilbert distance */
        int best_idx = 0;
        uint64_t best_diff = UINT64_MAX;
        for (uint32_t i = 0; i < node->count; i++) {
            uint64_t diff = (node->hilberts[i] > entry->hilbert) ?
                (node->hilberts[i] - entry->hilbert) :
                (entry->hilbert - node->hilberts[i]);
            if (diff < best_diff) {
                best_diff = diff;
                best_idx  = (int)i;
            }
        }

        /* If node is empty (fresh root), stop here */
        if (node->count == 0) break;

        node = node->children[best_idx];
    }

    return spatial_hrt_insert_at(rt, path, path_len,
                                  entry, entry->is_leaf, queue);
}

/* ============================================================================
 * Public API
 * ============================================================================ */

SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_new(void)
{
    return spatial_hilbertrtree_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();

    spatial_hilbertrtree *rt = (spatial_hilbertrtree *)alloc->malloc(
        sizeof(spatial_hilbertrtree), alloc->udata);
    if (SPATIAL_UNLIKELY(!rt)) return NULL;

    memset(rt, 0, sizeof(spatial_hilbertrtree));
    spatial_refcount_init(&rt->refcount);
    rt->alloc           = *alloc;
    rt->relaxed_atomics = true;

    rt->root = spatial_hrt_node_new(true, 0, alloc);
    if (SPATIAL_UNLIKELY(!rt->root)) {
        alloc->free(rt, alloc->udata);
        return NULL;
    }

    return rt;
}

SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_new_with_arena(
    spatial_arena_t *arena)
{
    spatial_allocator a = spatial_arena_allocator(arena);
    return spatial_hilbertrtree_new_with_allocator(&a);
}

SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_new_with_pool(
    spatial_pool_t *pool)
{
    spatial_allocator a = spatial_pool_allocator(pool);
    return spatial_hilbertrtree_new_with_allocator(&a);
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_new_with_memento_arena(
    memento_arena_t *arena)
{
    spatial_allocator a = spatial_memento_arena_allocator(arena);
    return spatial_hilbertrtree_new_with_allocator(&a);
}

SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_new_with_memento_pool(
    memento_pool_t *pool)
{
    spatial_allocator a = spatial_memento_pool_allocator(pool);
    return spatial_hilbertrtree_new_with_allocator(&a);
}
#endif

SPATIAL_INLINE void spatial_hilbertrtree_free(spatial_hilbertrtree *rt)
{
    if (SPATIAL_UNLIKELY(!rt)) return;
    if (spatial_refcount_decrement(&rt->refcount) > 0) return;

    spatial_hrt_node_free(rt->root, &rt->alloc, &rt->callbacks, rt->udata);
    rt->alloc.free(rt, rt->alloc.udata);
}

SPATIAL_INLINE spatial_hilbertrtree *spatial_hilbertrtree_clone(
    const spatial_hilbertrtree *rt)
{
    if (SPATIAL_UNLIKELY(!rt)) return NULL;
    spatial_refcount_increment(&rt->refcount);
    /* Cast away const: clone shares ownership, mutations must CoW */
    union { const spatial_hilbertrtree *c; spatial_hilbertrtree *m; } u;
    u.c = rt;
    return u.m;
}

SPATIAL_INLINE void spatial_hilbertrtree_set_item_callbacks(
    spatial_hilbertrtree *rt, const spatial_item_callbacks *cb)
{
    if (rt && cb) rt->callbacks = *cb;
}

SPATIAL_INLINE void spatial_hilbertrtree_set_udata(
    spatial_hilbertrtree *rt, void *udata)
{
    if (rt) rt->udata = udata;
}

SPATIAL_INLINE void spatial_hilbertrtree_opt_relaxed_atomics(
    spatial_hilbertrtree *rt, bool enable)
{
    if (rt) rt->relaxed_atomics = enable;
}

SPATIAL_INLINE bool spatial_hilbertrtree_insert(
    spatial_hilbertrtree *rt,
    const spatial_num_t *SPATIAL_RESTRICT min,
    const spatial_num_t *SPATIAL_RESTRICT max,
    spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!rt || !min || !max)) return false;

    /* Extend global bounds used for Hilbert normalization */
    if (!rt->bounds_valid) {
        spatial_bbox_copy(min, max, rt->global_min, rt->global_max,
                          SPATIAL_HILBERTRTREE_DIMS);
        rt->bounds_valid = true;
    } else {
        for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
            rt->global_min[d] = spatial_min(rt->global_min[d], min[d]);
            rt->global_max[d] = spatial_max(rt->global_max[d], max[d]);
        }
    }

    spatial_hrt_entry entry;
    memcpy(entry.min, min, sizeof(entry.min));
    memcpy(entry.max, max, sizeof(entry.max));
    entry.hilbert    = spatial_hrt_hilbert_for_tree(rt, min, max);
    entry.ptr.data   = data;
    entry.is_leaf    = true;

    spatial_hrt_reinsert_queue queue;
    memset(&queue, 0, sizeof(queue));

    bool ok = spatial_hrt_insert_entry(rt, &entry, &queue);

    /* Process any entries queued by the R* forced-reinsert step */
    for (int i = 0; i < queue.count && ok; i++) {
        queue.depth++;
        ok = spatial_hrt_insert_entry(rt, &queue.items[i].entry, &queue);
        queue.depth--;
    }

    if (ok) rt->count++;
    return ok;
}

SPATIAL_INLINE void spatial_hilbertrtree_search(
    const spatial_hilbertrtree *rt,
    const spatial_num_t *SPATIAL_RESTRICT qmin,
    const spatial_num_t *SPATIAL_RESTRICT qmax,
    spatial_hilbertrtree_iter_fn iter,
    void *udata)
{
    if (SPATIAL_UNLIKELY(!rt || !rt->root || !iter)) return;

    /* Fixed-size on-stack traversal — max height with branching 16 is tiny */
    spatial_hilbertrtree_node *stk[64];
    int sp = 0;
    stk[sp++] = rt->root;

    while (sp > 0) {
        spatial_hilbertrtree_node *node = stk[--sp];

        /* Node-level rejection */
        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(
                node->node_min, node->node_max, qmin, qmax,
                SPATIAL_HILBERTRTREE_DIMS))) {
            continue;
        }

        /* Prefetch first 8 children while we process entries */
        if (!node->is_leaf) {
            for (uint32_t i = 0; i < node->count && i < 8; i++) {
#if defined(__GNUC__) || defined(__clang__)
                __builtin_prefetch(node->children[i], 0, 1);
#endif
            }
        }

        for (uint32_t i = 0; i < node->count; i++) {
            if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(
                    node->mins[i], node->maxs[i], qmin, qmax,
                    SPATIAL_HILBERTRTREE_DIMS))) {
                continue;
            }

            if (node->is_leaf) {
                if (!iter(node->mins[i], node->maxs[i], node->datas[i], udata)) {
                    return;
                }
            } else {
                if (SPATIAL_LIKELY(sp < 64)) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE void spatial_hilbertrtree_scan(
    const spatial_hilbertrtree *rt,
    spatial_hilbertrtree_iter_fn iter,
    void *udata)
{
    if (SPATIAL_UNLIKELY(!rt || !rt->root || !iter)) return;

    spatial_hilbertrtree_node *stk[64];
    int sp = 0;
    stk[sp++] = rt->root;

    while (sp > 0) {
        spatial_hilbertrtree_node *node = stk[--sp];

        for (uint32_t i = 0; i < node->count; i++) {
            if (node->is_leaf) {
                if (!iter(node->mins[i], node->maxs[i], node->datas[i], udata)) {
                    return;
                }
            } else {
                if (SPATIAL_LIKELY(sp < 64)) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE int spatial_hilbertrtree_count(const spatial_hilbertrtree *rt)
{
    return rt ? rt->count : 0;
}

SPATIAL_INLINE bool spatial_hilbertrtree_delete(
    spatial_hilbertrtree *rt,
    const spatial_num_t *min,
    const spatial_num_t *max,
    spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!rt || !rt->root)) return false;

    spatial_hilbertrtree_node *stk[64];
    int sp = 0;
    stk[sp++] = rt->root;

    while (sp > 0) {
        spatial_hilbertrtree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->node_min, node->node_max, min, max,
                                   SPATIAL_HILBERTRTREE_DIMS)) {
            continue;
        }

        for (uint32_t i = 0; i < node->count; i++) {
            if (!spatial_bbox_overlaps(node->mins[i], node->maxs[i], min, max,
                                       SPATIAL_HILBERTRTREE_DIMS)) {
                continue;
            }

            if (node->is_leaf && node->datas[i] == data) {
                if (rt->callbacks.free) {
                    rt->callbacks.free(node->datas[i], rt->udata);
                }
                /* Compact by swapping with last entry */
                uint32_t last = node->count - 1;
                if (i != last) {
                    memcpy(node->mins[i],  node->mins[last],  sizeof(node->mins[0]));
                    memcpy(node->maxs[i],  node->maxs[last],  sizeof(node->maxs[0]));
                    node->hilberts[i]  = node->hilberts[last];
                    node->datas[i]     = node->datas[last];
                }
                node->count--;
                spatial_hrt_node_update_bounds(node);
                rt->count--;
                return true;
            } else if (!node->is_leaf) {
                if (SPATIAL_LIKELY(sp < 64)) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }

    return false;
}

/* ============================================================================
 * Bulk loading — sort by Hilbert value then STR-pack into the tree
 * ============================================================================ */

typedef struct {
    spatial_num_t min[SPATIAL_HILBERTRTREE_DIMS];
    spatial_num_t max[SPATIAL_HILBERTRTREE_DIMS];
    spatial_data_t data;
    uint64_t hilbert;
} spatial_hilbertrtree_bulk_item;

SPATIAL_INLINE int spatial_hrt_bulk_cmp(const void *a, const void *b)
{
    const spatial_hilbertrtree_bulk_item *ia = (const spatial_hilbertrtree_bulk_item *)a;
    const spatial_hilbertrtree_bulk_item *ib = (const spatial_hilbertrtree_bulk_item *)b;
    if (ia->hilbert < ib->hilbert) return -1;
    if (ia->hilbert > ib->hilbert) return  1;
    return 0;
}

SPATIAL_INLINE spatial_hilbertrtree_node *spatial_hrt_build_level(
    spatial_hilbertrtree_bulk_item *items, int count,
    int level, const spatial_hilbertrtree *rt);

SPATIAL_INLINE spatial_hilbertrtree_node *spatial_hrt_build_level(
    spatial_hilbertrtree_bulk_item *items, int count,
    int level, const spatial_hilbertrtree *rt)
{
    if (count <= 0) return NULL;

    spatial_hilbertrtree_node *node =
        spatial_hrt_node_new(level == 0, (uint32_t)level, &rt->alloc);
    if (SPATIAL_UNLIKELY(!node)) return NULL;

    if (level == 0) {
        int fill = count < SPATIAL_HILBERTRTREE_MAX_CHILDREN ?
                   count : SPATIAL_HILBERTRTREE_MAX_CHILDREN;
        for (int i = 0; i < fill; i++) {
            memcpy(node->mins[i],  items[i].min, sizeof(node->mins[i]));
            memcpy(node->maxs[i],  items[i].max, sizeof(node->maxs[i]));
            node->hilberts[i] = items[i].hilbert;
            node->datas[i]    = items[i].data;
        }
        node->count = (uint32_t)fill;
        spatial_hrt_node_update_bounds(node);
        return node;
    }

    /* Compute items_per_child = MAX_CHILDREN^(level) */
    int items_per_child = 1;
    for (int l = 0; l < level; l++) items_per_child *= SPATIAL_HILBERTRTREE_MAX_CHILDREN;

    int num_children = (count + items_per_child - 1) / items_per_child;
    if (num_children > SPATIAL_HILBERTRTREE_MAX_CHILDREN)
        num_children = SPATIAL_HILBERTRTREE_MAX_CHILDREN;

    for (int c = 0; c < num_children; c++) {
        int start = c * items_per_child;
        int end   = start + items_per_child;
        if (end > count) end = count;

        spatial_hilbertrtree_node *child =
            spatial_hrt_build_level(items + start, end - start, level - 1, rt);
        if (!child) continue;

        uint32_t idx = node->count;
        memcpy(node->mins[idx],  child->node_min, sizeof(node->mins[idx]));
        memcpy(node->maxs[idx],  child->node_max, sizeof(node->maxs[idx]));
        node->hilberts[idx]  = spatial_hrt_hilbert_for_tree(rt, child->node_min, child->node_max);
        node->children[idx]  = child;
        node->count++;
    }

    spatial_hrt_node_update_bounds(node);
    return node;
}

SPATIAL_INLINE bool spatial_hilbertrtree_bulk_load(
    spatial_hilbertrtree *rt,
    spatial_hilbertrtree_bulk_item *items,
    int count)
{
    if (SPATIAL_UNLIKELY(!rt || !items || count <= 0)) return false;

    /* Compute global bounds and encode Hilbert values */
    for (int i = 0; i < count; i++) {
        if (!rt->bounds_valid) {
            spatial_bbox_copy(items[i].min, items[i].max,
                              rt->global_min, rt->global_max,
                              SPATIAL_HILBERTRTREE_DIMS);
            rt->bounds_valid = true;
        } else {
            for (int d = 0; d < SPATIAL_HILBERTRTREE_DIMS; d++) {
                rt->global_min[d] = spatial_min(rt->global_min[d], items[i].min[d]);
                rt->global_max[d] = spatial_max(rt->global_max[d], items[i].max[d]);
            }
        }
    }
    for (int i = 0; i < count; i++) {
        items[i].hilbert = spatial_hrt_hilbert_for_tree(rt, items[i].min, items[i].max);
    }

    qsort(items, (size_t)count, sizeof(spatial_hilbertrtree_bulk_item),
          spatial_hrt_bulk_cmp);

    /* Compute required tree height */
    int cap = SPATIAL_HILBERTRTREE_MAX_CHILDREN;
    int height = 1;
    while (cap < count) {
        cap    *= SPATIAL_HILBERTRTREE_MAX_CHILDREN;
        height++;
    }

    spatial_hilbertrtree_node *new_root =
        spatial_hrt_build_level(items, count, height - 1, rt);

    if (SPATIAL_UNLIKELY(!new_root)) return false;

    spatial_hrt_node_free(rt->root, &rt->alloc, &rt->callbacks, rt->udata);
    rt->root   = new_root;
    rt->height = height;
    rt->count  = count;
    return true;
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_HILBERTRTREE_H */
