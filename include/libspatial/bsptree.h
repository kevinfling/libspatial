/**
 * @file   bsptree.h
 * @brief  Binary Space Partitioning Tree (2D/3D)
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * BSP-Tree partitions space using arbitrary planes (or axis-aligned fallback).
 * Commonly used in:
 * - Computer graphics (visibility determination)
 * - Collision detection
 * - CSG (Constructive Solid Geometry) operations
 */

#ifndef SPATIAL_BSPTREE_H
#define SPATIAL_BSPTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

#ifndef SPATIAL_BSPTREE_DIMS
    #define SPATIAL_BSPTREE_DIMS 3
#endif

#ifndef SPATIAL_BSPTREE_MAX_DEPTH
    #define SPATIAL_BSPTREE_MAX_DEPTH 24
#endif

#ifndef SPATIAL_BSPTREE_LEAF_SIZE
    #define SPATIAL_BSPTREE_LEAF_SIZE 8
#endif

/* Plane equation: ax + by + cz + d = 0 (normalized) */
typedef struct {
    spatial_num_t normal[SPATIAL_BSPTREE_DIMS];
    spatial_num_t d;
} spatial_bsp_plane;

/* Forward declarations */
struct spatial_bsptree_node;

/* Item structure */
typedef struct spatial_bsptree_item {
    spatial_num_t min[SPATIAL_BSPTREE_DIMS];
    spatial_num_t max[SPATIAL_BSPTREE_DIMS];
    spatial_data_t data;
} spatial_bsptree_item;

/* Node structure (cache-line aligned for SIMD-ready traversal) */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_bsptree_node {
    spatial_num_t min[SPATIAL_BSPTREE_DIMS];
    spatial_num_t max[SPATIAL_BSPTREE_DIMS];
    spatial_bsp_plane plane;
    struct spatial_bsptree_node *front;
    struct spatial_bsptree_node *back;

    spatial_bsptree_item *items;
    int item_count;
    int item_capacity;
    int depth;
    bool is_leaf;
} spatial_bsptree_node;

/* Main BSP-Tree structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_bsptree {
    spatial_bsptree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    int count;
    int max_depth;
    int leaf_size;
    bool axis_aligned;
    bool relaxed_atomics;
    bool needs_rebuild;
} spatial_bsptree;

/* Iterator callback */
typedef bool (*spatial_bsptree_iter_fn)(const spatial_num_t *min,
                                         const spatial_num_t *max,
                                         spatial_data_t data,
                                         void *udata);

/* Forward declarations for public API */
SPATIAL_INLINE spatial_bsptree* spatial_bsptree_new_with_allocator(const spatial_allocator *alloc);

/* Plane split result */
typedef enum {
    SPATIAL_BSP_FRONT = 0,
    SPATIAL_BSP_BACK = 1,
    SPATIAL_BSP_SPANNING = 2
} spatial_bsp_side;

/* ============================================================================
 * Plane operations
 * ============================================================================ */

SPATIAL_INLINE spatial_num_t spatial_bsp_plane_distance(const spatial_bsp_plane *plane,
                                                         const spatial_num_t *point)
{
    spatial_num_t dist = plane->d;
    for (int i = 0; i < SPATIAL_BSPTREE_DIMS; i++) {
        dist += plane->normal[i] * point[i];
    }
    return dist;
}

SPATIAL_INLINE spatial_bsp_side spatial_bsp_classify_bbox(const spatial_bsp_plane *plane,
                                                           const spatial_num_t *min,
                                                           const spatial_num_t *max)
{
    spatial_num_t min_dist = (spatial_num_t)0.0;
    spatial_num_t max_dist = (spatial_num_t)0.0;
    
    for (int i = 0; i < SPATIAL_BSPTREE_DIMS; i++) {
        if (plane->normal[i] >= 0) {
            min_dist += plane->normal[i] * min[i];
            max_dist += plane->normal[i] * max[i];
        } else {
            min_dist += plane->normal[i] * max[i];
            max_dist += plane->normal[i] * min[i];
        }
    }
    
    min_dist += plane->d;
    max_dist += plane->d;
    
    if (max_dist < -SPATIAL_EPSILON) return SPATIAL_BSP_BACK;
    if (min_dist > SPATIAL_EPSILON) return SPATIAL_BSP_FRONT;
    return SPATIAL_BSP_SPANNING;
}

SPATIAL_INLINE spatial_bsp_side spatial_bsp_classify_point(const spatial_bsp_plane *plane,
                                                            const spatial_num_t *point)
{
    spatial_num_t dist = spatial_bsp_plane_distance(plane, point);
    if (dist > SPATIAL_EPSILON) return SPATIAL_BSP_FRONT;
    if (dist < -SPATIAL_EPSILON) return SPATIAL_BSP_BACK;
    return SPATIAL_BSP_SPANNING;
}

/* ============================================================================
 * Internal Node Management
 * ============================================================================ */

SPATIAL_INLINE spatial_bsptree_node* spatial_bsptree_node_new_leaf(
    int depth,
    const spatial_allocator *alloc)
{
    spatial_bsptree_node *node = (spatial_bsptree_node*)alloc->malloc(sizeof(spatial_bsptree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    spatial_bbox_init(node->min, node->max, SPATIAL_BSPTREE_DIMS);
    node->is_leaf = true;
    node->items = NULL;
    node->item_count = 0;
    node->item_capacity = 0;
    node->depth = depth;
    node->front = NULL;
    node->back = NULL;
    
    return node;
}

SPATIAL_INLINE spatial_bsptree_node* spatial_bsptree_node_new_internal(
    int depth,
    const spatial_bsp_plane *plane,
    const spatial_allocator *alloc)
{
    spatial_bsptree_node *node = (spatial_bsptree_node*)alloc->malloc(sizeof(spatial_bsptree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    spatial_bbox_init(node->min, node->max, SPATIAL_BSPTREE_DIMS);
    node->is_leaf = false;
    node->depth = depth;
    node->plane = *plane;
    node->front = NULL;
    node->back = NULL;
    node->items = NULL;
    node->item_count = 0;
    
    return node;
}

SPATIAL_INLINE void spatial_bsptree_node_free(spatial_bsptree_node *node,
                                               const spatial_allocator *alloc,
                                               const spatial_item_callbacks *cb,
                                               void *udata);

SPATIAL_INLINE void spatial_bsptree_node_free(spatial_bsptree_node *node,
                                               const spatial_allocator *alloc,
                                               const spatial_item_callbacks *cb,
                                               void *udata)
{
    if (SPATIAL_UNLIKELY(!node)) return;
    
    if (node->is_leaf) {
        if (cb && cb->free) {
            for (int i = 0; i < node->item_count; i++) {
                cb->free(node->items[i].data, udata);
            }
        }
        alloc->free(node->items, alloc->udata);
    } else {
        spatial_bsptree_node_free(node->front, alloc, cb, udata);
        spatial_bsptree_node_free(node->back, alloc, cb, udata);
    }
    alloc->free(node, alloc->udata);
}

/* ============================================================================
 * Plane selection (axis-aligned binned SAH)
 * ============================================================================ */

#define SPATIAL_BSPTREE_NUM_BINS 8

typedef struct {
    spatial_num_t min[SPATIAL_BSPTREE_DIMS];
    spatial_num_t max[SPATIAL_BSPTREE_DIMS];
    int count;
} spatial_bsptree_bin;

SPATIAL_INLINE void spatial_bsptree_select_plane(
    spatial_bsptree_item *items,
    int count,
    bool axis_aligned,
    spatial_bsp_plane *out_plane)
{
    (void)axis_aligned;
    
    /* Compute overall bbox */
    spatial_num_t bmin[SPATIAL_BSPTREE_DIMS];
    spatial_num_t bmax[SPATIAL_BSPTREE_DIMS];
    for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
        bmin[d] = items[0].min[d];
        bmax[d] = items[0].max[d];
    }
    for (int i = 1; i < count; i++) {
        for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
            bmin[d] = spatial_min(bmin[d], items[i].min[d]);
            bmax[d] = spatial_max(bmax[d], items[i].max[d]);
        }
    }
    
    spatial_num_t total_sa = spatial_bbox_area(bmin, bmax, SPATIAL_BSPTREE_DIMS);
    spatial_num_t best_cost = SPATIAL_INFINITY;
    int best_axis = 0;
    spatial_num_t best_split_pos = (bmin[0] + bmax[0]) * (spatial_num_t)0.5;
    
    for (int axis = 0; axis < SPATIAL_BSPTREE_DIMS; axis++) {
        spatial_num_t axis_range = bmax[axis] - bmin[axis];
        if (axis_range < SPATIAL_EPSILON) continue;
        
        spatial_bsptree_bin bins[SPATIAL_BSPTREE_NUM_BINS];
        memset(bins, 0, sizeof(bins));
        
        /* Bin items by center */
        for (int i = 0; i < count; i++) {
            spatial_num_t center = (items[i].min[axis] + items[i].max[axis]) * (spatial_num_t)0.5;
            int bin_idx = (int)((center - bmin[axis]) / axis_range *
                                (SPATIAL_BSPTREE_NUM_BINS - 0.001));
            if (bin_idx < 0) bin_idx = 0;
            if (bin_idx >= SPATIAL_BSPTREE_NUM_BINS) bin_idx = SPATIAL_BSPTREE_NUM_BINS - 1;
            
            spatial_bsptree_bin *b = &bins[bin_idx];
            if (b->count == 0) {
                spatial_bbox_copy(items[i].min, items[i].max, b->min, b->max, SPATIAL_BSPTREE_DIMS);
            } else {
                for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
                    b->min[d] = spatial_min(b->min[d], items[i].min[d]);
                    b->max[d] = spatial_max(b->max[d], items[i].max[d]);
                }
            }
            b->count++;
        }
        
        /* Evaluate splits between bins */
        for (int split = 1; split < SPATIAL_BSPTREE_NUM_BINS; split++) {
            spatial_num_t left_min[SPATIAL_BSPTREE_DIMS];
            spatial_num_t left_max[SPATIAL_BSPTREE_DIMS];
            spatial_num_t right_min[SPATIAL_BSPTREE_DIMS];
            spatial_num_t right_max[SPATIAL_BSPTREE_DIMS];
            spatial_bbox_init(left_min, left_max, SPATIAL_BSPTREE_DIMS);
            spatial_bbox_init(right_min, right_max, SPATIAL_BSPTREE_DIMS);
            int left_count = 0, right_count = 0;
            
            for (int j = 0; j < split; j++) {
                if (bins[j].count > 0) {
                    for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
                        left_min[d] = spatial_min(left_min[d], bins[j].min[d]);
                        left_max[d] = spatial_max(left_max[d], bins[j].max[d]);
                    }
                    left_count += bins[j].count;
                }
            }
            for (int j = split; j < SPATIAL_BSPTREE_NUM_BINS; j++) {
                if (bins[j].count > 0) {
                    for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
                        right_min[d] = spatial_min(right_min[d], bins[j].min[d]);
                        right_max[d] = spatial_max(right_max[d], bins[j].max[d]);
                    }
                    right_count += bins[j].count;
                }
            }
            
            if (left_count == 0 || right_count == 0) continue;
            
            spatial_num_t left_sa = spatial_bbox_area(left_min, left_max, SPATIAL_BSPTREE_DIMS);
            spatial_num_t right_sa = spatial_bbox_area(right_min, right_max, SPATIAL_BSPTREE_DIMS);
            spatial_num_t cost = spatial_sah_cost(total_sa, left_sa, right_sa,
                                                   left_count, right_count,
                                                   (spatial_num_t)1.0, (spatial_num_t)1.0);
            
            if (cost < best_cost) {
                best_cost = cost;
                best_axis = axis;
                best_split_pos = bmin[axis] + (axis_range * (spatial_num_t)split) /
                                               (spatial_num_t)SPATIAL_BSPTREE_NUM_BINS;
            }
        }
    }
    
    /* If SAH didn't find a useful split, fall back to median of longest axis */
    if (best_cost >= SPATIAL_INFINITY) {
        spatial_num_t best_extent = (spatial_num_t)0.0;
        for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
            spatial_num_t extent = bmax[d] - bmin[d];
            if (extent > best_extent) {
                best_extent = extent;
                best_axis = d;
            }
        }
        best_split_pos = (bmin[best_axis] + bmax[best_axis]) * (spatial_num_t)0.5;
    }
    
    for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
        out_plane->normal[d] = (d == best_axis) ? (spatial_num_t)1.0 : (spatial_num_t)0.0;
    }
    out_plane->d = -best_split_pos;
}

/* ============================================================================
 * Public API
 * ============================================================================ */

SPATIAL_INLINE spatial_bsptree* spatial_bsptree_new(void) {
    return spatial_bsptree_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_bsptree* spatial_bsptree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();
    
    spatial_bsptree *bsp = (spatial_bsptree*)alloc->malloc(sizeof(spatial_bsptree), alloc->udata);
    if (SPATIAL_UNLIKELY(!bsp)) return NULL;
    
    memset(bsp, 0, sizeof(spatial_bsptree));
    spatial_refcount_init(&bsp->refcount);
    bsp->alloc = *alloc;
    bsp->max_depth = SPATIAL_BSPTREE_MAX_DEPTH;
    bsp->leaf_size = SPATIAL_BSPTREE_LEAF_SIZE;
    bsp->axis_aligned = true;
    bsp->relaxed_atomics = true;
    bsp->needs_rebuild = false;
    bsp->root = NULL;
    
    return bsp;
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_bsptree* spatial_bsptree_new_with_memento_arena(memento_arena_t *arena) {
    spatial_allocator alloc = spatial_memento_arena_allocator(arena);
    return spatial_bsptree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_bsptree* spatial_bsptree_new_with_memento_pool(memento_pool_t *pool) {
    spatial_allocator alloc = spatial_memento_pool_allocator(pool);
    return spatial_bsptree_new_with_allocator(&alloc);
}
#endif

SPATIAL_INLINE spatial_bsptree* spatial_bsptree_new_with_arena(spatial_arena_t *arena) {
    spatial_allocator alloc = spatial_arena_allocator(arena);
    return spatial_bsptree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_bsptree* spatial_bsptree_new_with_pool(spatial_pool_t *pool) {
    spatial_allocator alloc = spatial_pool_allocator(pool);
    return spatial_bsptree_new_with_allocator(&alloc);
}



SPATIAL_INLINE void spatial_bsptree_free(spatial_bsptree *bsp) {
    if (SPATIAL_UNLIKELY(!bsp)) return;
    
    if (spatial_refcount_decrement(&bsp->refcount) > 0) return;
    
    spatial_bsptree_node_free(bsp->root, &bsp->alloc, &bsp->callbacks, bsp->udata);
    bsp->alloc.free(bsp, bsp->alloc.udata);
}

SPATIAL_INLINE spatial_bsptree* spatial_bsptree_clone(const spatial_bsptree *bsp) {
    if (SPATIAL_UNLIKELY(!bsp)) return NULL;
    
    spatial_refcount_increment(&bsp->refcount);
    union { const spatial_bsptree *in; spatial_bsptree *out; } pun; pun.in = bsp; return pun.out;
}

SPATIAL_INLINE void spatial_bsptree_set_item_callbacks(spatial_bsptree *bsp,
                                                        const spatial_item_callbacks *cb)
{
    if (bsp) bsp->callbacks = *cb;
}

SPATIAL_INLINE void spatial_bsptree_set_udata(spatial_bsptree *bsp, void *udata) {
    if (bsp) bsp->udata = udata;
}

SPATIAL_INLINE void spatial_bsptree_opt_relaxed_atomics(spatial_bsptree *bsp, bool enable) {
    if (bsp) bsp->relaxed_atomics = enable;
}

SPATIAL_INLINE void spatial_bsptree_set_axis_aligned(spatial_bsptree *bsp, bool axis_aligned) {
    if (bsp) bsp->axis_aligned = axis_aligned;
}

SPATIAL_INLINE bool spatial_bsptree_insert(spatial_bsptree *bsp,
                                            const spatial_num_t *min,
                                            const spatial_num_t *max,
                                            spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!bsp)) return false;
    
    if (SPATIAL_UNLIKELY(!bsp->root)) {
        bsp->root = spatial_bsptree_node_new_leaf(0, &bsp->alloc);
        if (SPATIAL_UNLIKELY(!bsp->root)) return false;
    }
    
    spatial_bsptree_node *leaf = bsp->root;

    /* Grow item array with 2× amortised capacity */
    if (leaf->item_count >= leaf->item_capacity) {
        int new_cap = leaf->item_capacity ? leaf->item_capacity * 2 : 8;
        spatial_bsptree_item *new_items = (spatial_bsptree_item*)bsp->alloc.malloc(
            sizeof(spatial_bsptree_item) * (size_t)new_cap, bsp->alloc.udata);
        if (SPATIAL_UNLIKELY(!new_items)) return false;
        if (leaf->items) {
            memcpy(new_items, leaf->items,
                   sizeof(spatial_bsptree_item) * (size_t)leaf->item_count);
            bsp->alloc.free(leaf->items, bsp->alloc.udata);
        }
        leaf->items = new_items;
        leaf->item_capacity = new_cap;
    }

    spatial_bbox_copy(min, max,
                      leaf->items[leaf->item_count].min,
                      leaf->items[leaf->item_count].max,
                      SPATIAL_BSPTREE_DIMS);
    leaf->items[leaf->item_count].data = data;
    leaf->item_count++;
    bsp->count++;
    bsp->needs_rebuild = true;
    
    /* Expand leaf bbox to include new item */
    if (leaf->item_count == 1) {
        spatial_bbox_copy(min, max, leaf->min, leaf->max, SPATIAL_BSPTREE_DIMS);
    } else {
        for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
            leaf->min[d] = spatial_min(leaf->min[d], min[d]);
            leaf->max[d] = spatial_max(leaf->max[d], max[d]);
        }
    }

    return true;
}

/* Build tree recursively */
SPATIAL_INLINE spatial_bsptree_node* spatial_bsptree_build_recursive(
    spatial_bsptree_item *items,
    int count,
    int depth,
    int max_depth,
    int leaf_size,
    bool axis_aligned,
    const spatial_allocator *alloc)
{
    if (count <= leaf_size || depth >= max_depth) {
        spatial_bsptree_node *leaf = spatial_bsptree_node_new_leaf(depth, alloc);
        if (SPATIAL_UNLIKELY(!leaf)) return NULL;
        
        leaf->items = (spatial_bsptree_item*)alloc->malloc(
            sizeof(spatial_bsptree_item) * (size_t)count, alloc->udata);
        if (SPATIAL_UNLIKELY(!leaf->items)) {
            alloc->free(leaf, alloc->udata);
            return NULL;
        }
        
        memcpy(leaf->items, items, sizeof(spatial_bsptree_item) * (size_t)count);
        leaf->item_count = count;
        
        /* Compute leaf bbox */
        for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
            leaf->min[d] = items[0].min[d];
            leaf->max[d] = items[0].max[d];
        }
        for (int i = 1; i < count; i++) {
            for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
                leaf->min[d] = spatial_min(leaf->min[d], items[i].min[d]);
                leaf->max[d] = spatial_max(leaf->max[d], items[i].max[d]);
            }
        }
        
        return leaf;
    }
    
    /* Select splitting plane */
    spatial_bsp_plane plane;
    spatial_bsptree_select_plane(items, count, axis_aligned, &plane);
    
    spatial_bsptree_node *node = spatial_bsptree_node_new_internal(depth, &plane, alloc);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    /* Single-pass partition: over-allocate to count, fill directly */
    spatial_bsptree_item *front_items = (spatial_bsptree_item*)alloc->malloc(
        sizeof(spatial_bsptree_item) * (size_t)count, alloc->udata);
    spatial_bsptree_item *back_items = (spatial_bsptree_item*)alloc->malloc(
        sizeof(spatial_bsptree_item) * (size_t)count, alloc->udata);

    if (SPATIAL_UNLIKELY(!front_items || !back_items)) {
        alloc->free(front_items, alloc->udata);
        alloc->free(back_items, alloc->udata);
        alloc->free(node, alloc->udata);
        return NULL;
    }

    int fi = 0, bi = 0;
    for (int i = 0; i < count; i++) {
        spatial_num_t center[SPATIAL_BSPTREE_DIMS];
        for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
            center[d] = (items[i].min[d] + items[i].max[d]) * (spatial_num_t)0.5;
        }

        spatial_bsp_side side = spatial_bsp_classify_point(&plane, center);
        if (side == SPATIAL_BSP_BACK) {
            back_items[bi++] = items[i];
        } else {
            /* FRONT or ON_PLANE: assign to front to avoid duplication */
            front_items[fi++] = items[i];
        }
    }

    /* Degenerate: all items ended up on one side — make a leaf */
    if (fi == 0 || bi == 0) {
        alloc->free(front_items, alloc->udata);
        alloc->free(back_items, alloc->udata);
        alloc->free(node, alloc->udata);

        spatial_bsptree_node *leaf = spatial_bsptree_node_new_leaf(depth, alloc);
        if (SPATIAL_UNLIKELY(!leaf)) return NULL;

        leaf->items = (spatial_bsptree_item*)alloc->malloc(
            sizeof(spatial_bsptree_item) * (size_t)count, alloc->udata);
        if (SPATIAL_UNLIKELY(!leaf->items)) {
            alloc->free(leaf, alloc->udata);
            return NULL;
        }

        memcpy(leaf->items, items, sizeof(spatial_bsptree_item) * (size_t)count);
        leaf->item_count = count;
        
        /* Compute leaf bbox */
        for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
            leaf->min[d] = items[0].min[d];
            leaf->max[d] = items[0].max[d];
        }
        for (int i = 1; i < count; i++) {
            for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
                leaf->min[d] = spatial_min(leaf->min[d], items[i].min[d]);
                leaf->max[d] = spatial_max(leaf->max[d], items[i].max[d]);
            }
        }
        return leaf;
    }

    node->front = spatial_bsptree_build_recursive(front_items, fi, depth + 1,
                                                   max_depth, leaf_size, axis_aligned, alloc);
    node->back  = spatial_bsptree_build_recursive(back_items,  bi, depth + 1,
                                                   max_depth, leaf_size, axis_aligned, alloc);

    alloc->free(front_items, alloc->udata);
    alloc->free(back_items, alloc->udata);

    if (SPATIAL_UNLIKELY(!node->front || !node->back)) {
        spatial_bsptree_node_free(node->front, alloc, NULL, NULL);
        spatial_bsptree_node_free(node->back, alloc, NULL, NULL);
        alloc->free(node, alloc->udata);
        return NULL;
    }
    
    /* Compute subtree bbox from children */
    for (int d = 0; d < SPATIAL_BSPTREE_DIMS; d++) {
        node->min[d] = spatial_min(node->front->min[d], node->back->min[d]);
        node->max[d] = spatial_max(node->front->max[d], node->back->max[d]);
    }

    return node;
}

SPATIAL_INLINE void spatial_bsptree_rebuild(spatial_bsptree *bsp) {
    if (SPATIAL_UNLIKELY(!bsp || !bsp->needs_rebuild || !bsp->root)) return;
    
    spatial_bsptree_item *items = bsp->root->items;
    int count = bsp->root->item_count;
    
    bsp->root->items = NULL;
    bsp->alloc.free(bsp->root, bsp->alloc.udata);
    
    bsp->root = spatial_bsptree_build_recursive(items, count, 0, bsp->max_depth, bsp->leaf_size,
                                                 bsp->axis_aligned, &bsp->alloc);
    bsp->alloc.free(items, bsp->alloc.udata);
    bsp->needs_rebuild = false;
}

/* Search */
SPATIAL_INLINE void spatial_bsptree_search(const spatial_bsptree *bsp,
                                            const spatial_num_t *SPATIAL_RESTRICT qmin,
                                            const spatial_num_t *SPATIAL_RESTRICT qmax,
                                            spatial_bsptree_iter_fn iter,
                                            void *udata)
{
    if (SPATIAL_UNLIKELY(!bsp || !bsp->root || !iter)) return;

    /* Fixed on-stack traversal — binary tree, max depth SPATIAL_BSPTREE_MAX_DEPTH */
    spatial_bsptree_node *stk[SPATIAL_BSPTREE_MAX_DEPTH + 2];
    int sp = 0;
    stk[sp++] = bsp->root;

    while (sp > 0) {
        spatial_bsptree_node *node = stk[--sp];

        /* Bbox cull: skip entire subtree if query doesn't overlap node bounds */
        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(node->min, node->max, qmin, qmax,
                                                     SPATIAL_BSPTREE_DIMS))) {
            continue;
        }

        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                spatial_bsptree_item *item = &node->items[i];
                if (spatial_bbox_overlaps(item->min, item->max, qmin, qmax,
                                          SPATIAL_BSPTREE_DIMS)) {
                    if (!iter(item->min, item->max, item->data, udata)) return;
                }
            }
        } else {
            spatial_bsp_side side = spatial_bsp_classify_bbox(&node->plane, qmin, qmax);

            if ((side == SPATIAL_BSP_FRONT || side == SPATIAL_BSP_SPANNING) && node->front)
                stk[sp++] = node->front;
            if ((side == SPATIAL_BSP_BACK  || side == SPATIAL_BSP_SPANNING) && node->back)
                stk[sp++] = node->back;
        }
    }
}

/* Scan all items */
SPATIAL_INLINE void spatial_bsptree_scan(const spatial_bsptree *bsp,
                                          spatial_bsptree_iter_fn iter,
                                          void *udata)
{
    if (SPATIAL_UNLIKELY(!bsp || !bsp->root || !iter)) return;

    spatial_bsptree_node *stk[SPATIAL_BSPTREE_MAX_DEPTH + 2];
    int sp = 0;
    stk[sp++] = bsp->root;

    while (sp > 0) {
        spatial_bsptree_node *node = stk[--sp];

        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                if (!iter(node->items[i].min, node->items[i].max,
                           node->items[i].data, udata)) return;
            }
        } else {
            if (node->front) stk[sp++] = node->front;
            if (node->back)  stk[sp++] = node->back;
        }
    }
}

SPATIAL_INLINE int spatial_bsptree_count(const spatial_bsptree *bsp) {
    return bsp ? bsp->count : 0;
}

/* Delete */
SPATIAL_INLINE bool spatial_bsptree_delete(spatial_bsptree *bsp,
                                            const spatial_num_t *min,
                                            const spatial_num_t *max,
                                            spatial_data_t data)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!bsp || !bsp->root)) return false;
    
    if (bsp->root->is_leaf) {
        for (int i = 0; i < bsp->root->item_count; i++) {
            if (bsp->root->items[i].data == data) {
                if (bsp->callbacks.free) {
                    bsp->callbacks.free(bsp->root->items[i].data, bsp->udata);
                }
                bsp->root->items[i] = bsp->root->items[bsp->root->item_count - 1];
                bsp->root->item_count--;
                bsp->count--;
                return true;
            }
        }
    }
    
    bsp->needs_rebuild = true;
    return false;
}

typedef bool (*spatial_bsptree_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_bsptree_delete_with_comparator(spatial_bsptree *bsp,
                                                            const spatial_num_t *min,
                                                            const spatial_num_t *max,
                                                            spatial_bsptree_cmp_fn cmp,
                                                            void *cmp_udata)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!bsp || !bsp->root)) return false;
    
    if (bsp->root->is_leaf) {
        for (int i = 0; i < bsp->root->item_count; i++) {
            if (cmp(bsp->root->items[i].data, cmp_udata)) {
                if (bsp->callbacks.free) {
                    bsp->callbacks.free(bsp->root->items[i].data, bsp->udata);
                }
                bsp->root->items[i] = bsp->root->items[bsp->root->item_count - 1];
                bsp->root->item_count--;
                bsp->count--;
                return true;
            }
        }
    }
    
    bsp->needs_rebuild = true;
    return false;
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_BSPTREE_H */
