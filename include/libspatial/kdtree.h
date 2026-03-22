/**
 * @file   kdtree.h
 * @brief  N-Dimensional KD-Tree with incremental splits
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * Incremental insertion with on-the-fly leaf splitting.
 * No lazy rebuilds. Median/sliding-midpoint for hot path,
 * SAH available for explicit static rebuilds.
 */

#ifndef SPATIAL_KDTREE_H
#define SPATIAL_KDTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

#ifndef SPATIAL_KDTREE_DIMS
    #define SPATIAL_KDTREE_DIMS 3
#endif

#ifndef SPATIAL_KDTREE_SPLIT
    #define SPATIAL_KDTREE_SPLIT SPATIAL_KDTREE_SPLIT_MEDIAN
#endif

#define SPATIAL_KDTREE_SPLIT_MEDIAN       0
#define SPATIAL_KDTREE_SPLIT_SAH          1
#define SPATIAL_KDTREE_SPLIT_SLIDING_MID  2

#ifndef SPATIAL_KDTREE_MAX_DEPTH
    #define SPATIAL_KDTREE_MAX_DEPTH 32
#endif

#ifndef SPATIAL_KDTREE_LEAF_SIZE
    #define SPATIAL_KDTREE_LEAF_SIZE 8
#endif

_Static_assert(SPATIAL_KDTREE_DIMS > 0, "KD-Tree requires positive dimensionality");
_Static_assert(SPATIAL_KDTREE_LEAF_SIZE >= 1, "Leaf size must be at least 1");

/* Forward declarations */
struct spatial_kdtree_node;

/* Item (point) */
typedef struct spatial_kdtree_item {
    spatial_num_t point[SPATIAL_KDTREE_DIMS];
    spatial_data_t data;
} spatial_kdtree_item;

/* Node — cache-line aligned, restrict-ready */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_kdtree_node {
    spatial_num_t min[SPATIAL_KDTREE_DIMS];
    spatial_num_t max[SPATIAL_KDTREE_DIMS];

    union {
        struct {
            struct spatial_kdtree_node *SPATIAL_RESTRICT left;
            struct spatial_kdtree_node *SPATIAL_RESTRICT right;
            spatial_num_t split_pos;
            int split_axis;
        };
        struct {
            spatial_kdtree_item *SPATIAL_RESTRICT items;
            int item_count;
            int item_capacity;
        };
    };

    bool is_leaf;
} spatial_kdtree_node;

/* Main structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_kdtree {
    spatial_kdtree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    int count;
    int max_depth;
    int leaf_size;
    bool relaxed_atomics;
} spatial_kdtree;

/* Iterator / distance callbacks */
typedef bool (*spatial_kdtree_iter_fn)(const spatial_num_t *min,
                                        const spatial_num_t *max,
                                        spatial_data_t data,
                                        void *udata);

typedef spatial_num_t (*spatial_kdtree_dist_fn)(const spatial_num_t *a,
                                                 const spatial_num_t *b,
                                                 int dims,
                                                 void *udata);

/* Forward declarations */
SPATIAL_INLINE spatial_kdtree* spatial_kdtree_new_with_allocator(const spatial_allocator *alloc);

SPATIAL_INLINE spatial_kdtree* spatial_kdtree_new(void) {
    return spatial_kdtree_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_kdtree* spatial_kdtree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();

    spatial_kdtree *kt = (spatial_kdtree*)alloc->malloc(sizeof(*kt), alloc->udata);
    if (SPATIAL_UNLIKELY(!kt)) return NULL;
    memset(kt, 0, sizeof(*kt));
    spatial_refcount_init(&kt->refcount);
    kt->alloc = *alloc;
    kt->max_depth = SPATIAL_KDTREE_MAX_DEPTH;
    kt->leaf_size = SPATIAL_KDTREE_LEAF_SIZE;
    kt->relaxed_atomics = true;
    return kt;
}

/* Default distance functions */
SPATIAL_INLINE SPATIAL_CONST spatial_num_t spatial_kdtree_dist_sq_euclidean(
    const spatial_num_t *a, const spatial_num_t *b, int dims, void *udata)
{
    (void)udata;
    spatial_num_t sum = (spatial_num_t)0.0;
    for (int i = 0; i < dims; i++) {
        spatial_num_t diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum;
}

SPATIAL_INLINE SPATIAL_CONST spatial_num_t spatial_kdtree_dist_manhattan(
    const spatial_num_t *a, const spatial_num_t *b, int dims, void *udata)
{
    (void)udata;
    spatial_num_t sum = (spatial_num_t)0.0;
    for (int i = 0; i < dims; i++) sum += spatial_abs(a[i] - b[i]);
    return sum;
}

/* Node management (memento-aware) */
SPATIAL_INLINE spatial_kdtree_node* spatial_kdtree_node_new_leaf(const spatial_allocator *alloc)
{
    spatial_kdtree_node *node = (spatial_kdtree_node*)alloc->malloc(sizeof(*node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    memset(node, 0, sizeof(*node));
    node->is_leaf = true;
    return node;
}

SPATIAL_INLINE spatial_kdtree_node* spatial_kdtree_node_new_internal(
    int split_axis, spatial_num_t split_pos, const spatial_allocator *alloc)
{
    spatial_kdtree_node *node = (spatial_kdtree_node*)alloc->malloc(sizeof(*node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    memset(node, 0, sizeof(*node));
    node->is_leaf = false;
    node->split_axis = split_axis;
    node->split_pos = split_pos;
    return node;
}

SPATIAL_INLINE void spatial_kdtree_node_free(spatial_kdtree_node *node,
                                              const spatial_allocator *alloc,
                                              const spatial_item_callbacks *cb,
                                              void *udata)
{
    if (SPATIAL_UNLIKELY(!node)) return;

    spatial_kdtree_node *stk[SPATIAL_KDTREE_MAX_DEPTH * 2];
    int sp = 0;
    stk[sp++] = node;

    while (sp > 0) {
        spatial_kdtree_node *n = stk[--sp];

        if (n->is_leaf) {
            if (cb && cb->free) {
                for (int i = 0; i < n->item_count; i++) {
                    cb->free(n->items[i].data, udata);
                }
            }
            alloc->free(n->items, alloc->udata);
        } else {
            if (n->left)  stk[sp++] = n->left;
            if (n->right) stk[sp++] = n->right;
        }
        alloc->free(n, alloc->udata);
    }
}

/* Incremental split helpers */
SPATIAL_INLINE int spatial_kdtree_choose_split_axis(const spatial_num_t *min,
                                                     const spatial_num_t *max,
                                                     int dims)
{
    int axis = 0;
    spatial_num_t max_extent = max[0] - min[0];
    for (int i = 1; i < dims; ++i) {
        spatial_num_t e = max[i] - min[i];
        if (e > max_extent) { max_extent = e; axis = i; }
    }
    return axis;
}

SPATIAL_INLINE spatial_num_t spatial_kdtree_split_median(spatial_kdtree_item *items,
                                                         int count, int axis)
{
    for (int i = 1; i < count; ++i) {
        spatial_kdtree_item key = items[i];
        int j = i - 1;
        while (j >= 0 && items[j].point[axis] > key.point[axis]) {
            items[j + 1] = items[j]; --j;
        }
        items[j + 1] = key;
    }
    return items[count / 2].point[axis];
}

SPATIAL_INLINE spatial_num_t spatial_kdtree_split_sliding_mid(const spatial_num_t *min,
                                                               const spatial_num_t *max, int axis)
{
    return (min[axis] + max[axis]) * (spatial_num_t)0.5;
}

SPATIAL_INLINE bool spatial_kdtree_node_split(spatial_kdtree_node *node,
                                              int depth,
                                              const spatial_allocator *alloc)
{
    if (depth >= SPATIAL_KDTREE_MAX_DEPTH) return false;

    int axis = spatial_kdtree_choose_split_axis(node->min, node->max, SPATIAL_KDTREE_DIMS);
    spatial_num_t split_pos = (SPATIAL_KDTREE_SPLIT == SPATIAL_KDTREE_SPLIT_SLIDING_MID)
        ? spatial_kdtree_split_sliding_mid(node->min, node->max, axis)
        : spatial_kdtree_split_median(node->items, node->item_count, axis);

    spatial_kdtree_node *left  = spatial_kdtree_node_new_internal(axis, split_pos, alloc);
    spatial_kdtree_node *right = spatial_kdtree_node_new_internal(axis, split_pos, alloc);
    if (SPATIAL_UNLIKELY(!left || !right)) return false;

    int left_count = 0;
    for (int i = 0; i < node->item_count; ++i)
        if (node->items[i].point[axis] <= split_pos) ++left_count;

    spatial_kdtree_item *left_items  = (spatial_kdtree_item*)alloc->malloc(
        sizeof(spatial_kdtree_item) * (size_t)left_count, alloc->udata);
    spatial_kdtree_item *right_items = (spatial_kdtree_item*)alloc->malloc(
        sizeof(spatial_kdtree_item) * (size_t)(node->item_count - left_count), alloc->udata);

    if (SPATIAL_UNLIKELY(!left_items || !right_items)) return false;

    int li = 0, ri = 0;
    for (int i = 0; i < node->item_count; ++i) {
        if (node->items[i].point[axis] <= split_pos)
            left_items[li++] = node->items[i];
        else
            right_items[ri++] = node->items[i];
    }

    left->is_leaf = true;  left->items = left_items;  left->item_count = left_count;  left->item_capacity = left_count;
    right->is_leaf = true; right->items = right_items; right->item_count = node->item_count - left_count; right->item_capacity = node->item_count - left_count;

    spatial_bbox_copy(node->min, node->max, left->min, left->max, SPATIAL_KDTREE_DIMS);
    spatial_bbox_copy(node->min, node->max, right->min, right->max, SPATIAL_KDTREE_DIMS);
    left->max[axis] = split_pos;
    right->min[axis] = split_pos;

    /* --- FIX: Cache items pointer before overwriting union --- */
    spatial_kdtree_item *old_items = node->items;

    node->is_leaf = false;
    node->left = left;
    node->right = right;
    node->split_axis = axis;
    node->split_pos = split_pos;

    alloc->free(old_items, alloc->udata);
    return true;
}

/* Incremental insert — the hot path */
SPATIAL_INLINE bool spatial_kdtree_insert(spatial_kdtree *kt,
                                           const spatial_num_t *SPATIAL_RESTRICT min,
                                           const spatial_num_t *SPATIAL_RESTRICT max,
                                           spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!kt)) return false;

    spatial_num_t point[SPATIAL_KDTREE_DIMS];
    for (int i = 0; i < SPATIAL_KDTREE_DIMS; ++i)
        point[i] = (min[i] + max[i]) * (spatial_num_t)0.5;

    if (SPATIAL_UNLIKELY(!kt->root)) {
        kt->root = spatial_kdtree_node_new_leaf(&kt->alloc);
        if (SPATIAL_UNLIKELY(!kt->root)) return false;
        spatial_bbox_init(kt->root->min, kt->root->max, SPATIAL_KDTREE_DIMS);
    }

    /* Record path for upward bounds propagation */
    spatial_kdtree_node *path[SPATIAL_KDTREE_MAX_DEPTH + 2];
    int path_len = 0;

    spatial_kdtree_node *node = kt->root;
    int depth = 0;

    while (!node->is_leaf) {
        path[path_len++] = node;
        int axis = node->split_axis;
        node = (point[axis] <= node->split_pos) ? node->left : node->right;
        ++depth;
    }

    /* Add to leaf with exponential capacity growth */
    int new_count = node->item_count + 1;
    if (new_count > node->item_capacity) {
        int new_cap = node->item_capacity ? node->item_capacity * 2 : kt->leaf_size;
        if (new_cap < new_count) new_cap = new_count;
        spatial_kdtree_item *new_items = (spatial_kdtree_item*)kt->alloc.malloc(
            sizeof(spatial_kdtree_item) * (size_t)new_cap, kt->alloc.udata);
        if (SPATIAL_UNLIKELY(!new_items)) return false;

        if (node->items) {
            memcpy(new_items, node->items, sizeof(spatial_kdtree_item) * (size_t)node->item_count);
            kt->alloc.free(node->items, kt->alloc.udata);
        }
        node->items = new_items;
        node->item_capacity = new_cap;
    }

    for (int i = 0; i < SPATIAL_KDTREE_DIMS; ++i) {
        node->items[node->item_count].point[i] = point[i];
        /* Expand leaf bounding box to include new point */
        node->min[i] = spatial_min(node->min[i], point[i]);
        node->max[i] = spatial_max(node->max[i], point[i]);
    }
    node->items[node->item_count].data = data;

    node->item_count = new_count;
    kt->count++;

    /* Split if needed */
    if (node->item_count > kt->leaf_size && depth < kt->max_depth) {
        spatial_kdtree_node_split(node, depth, &kt->alloc);
    }

    /* Propagate bounds upward so internal node bounds stay accurate */
    for (int pi = path_len - 1; pi >= 0; pi--) {
        spatial_kdtree_node *p = path[pi];
        bool changed = false;
        for (int i = 0; i < SPATIAL_KDTREE_DIMS; ++i) {
            if (point[i] < p->min[i]) { p->min[i] = point[i]; changed = true; }
            if (point[i] > p->max[i]) { p->max[i] = point[i]; changed = true; }
        }
        if (!changed) break; /* no further ancestors need updating */
    }

    return true;
}

/* Search */
SPATIAL_INLINE void spatial_kdtree_search(const spatial_kdtree *kt,
                                           const spatial_num_t *SPATIAL_RESTRICT qmin,
                                           const spatial_num_t *SPATIAL_RESTRICT qmax,
                                           spatial_kdtree_iter_fn iter,
                                           void *udata)
{
    if (SPATIAL_UNLIKELY(!kt || !kt->root || !iter)) return;

    spatial_kdtree_node *stk[SPATIAL_KDTREE_MAX_DEPTH + 2];
    int sp = 0;
    stk[sp++] = kt->root;

    while (sp > 0) {
        spatial_kdtree_node *node = stk[--sp];

        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(node->min, node->max, qmin, qmax,
                                                     SPATIAL_KDTREE_DIMS))) {
            continue;
        }

        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                spatial_kdtree_item *item = &node->items[i];
                if (spatial_bbox_contains(qmin, qmax, item->point, item->point,
                                          SPATIAL_KDTREE_DIMS)) {
                    if (!iter(item->point, item->point, item->data, udata)) return;
                }
            }
        } else {
            int axis = node->split_axis;
            spatial_num_t pos = node->split_pos;
            if (qmin[axis] <= pos && node->left) {
                SPATIAL_PREFETCH(node->left, 0, 3);
                stk[sp++] = node->left;
            }
            if (qmax[axis] >= pos && node->right) {
                SPATIAL_PREFETCH(node->right, 0, 3);
                stk[sp++] = node->right;
            }
        }
    }
}

/* Scan all items */
SPATIAL_INLINE void spatial_kdtree_scan(const spatial_kdtree *kt,
                                         spatial_kdtree_iter_fn iter,
                                         void *udata)
{
    if (SPATIAL_UNLIKELY(!kt || !kt->root || !iter)) return;

    spatial_kdtree_node *stk[SPATIAL_KDTREE_MAX_DEPTH + 2];
    int sp = 0;
    stk[sp++] = kt->root;

    while (sp > 0) {
        spatial_kdtree_node *node = stk[--sp];

        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                spatial_kdtree_item *item = &node->items[i];
                if (!iter(item->point, item->point, item->data, udata)) return;
            }
        } else {
            if (node->left)  stk[sp++] = node->left;
            if (node->right) stk[sp++] = node->right;
        }
    }
}

SPATIAL_INLINE int spatial_kdtree_count(const spatial_kdtree *kt) {
    return kt ? kt->count : 0;
}

/* k-NN search */
SPATIAL_INLINE void spatial_kdtree_nearest(
    const spatial_kdtree *kt,
    const spatial_num_t *query,
    int k,
    spatial_kdtree_dist_fn dist_fn,
    void *dist_udata,
    void **results,
    spatial_num_t *distances)
{
    if (SPATIAL_UNLIKELY(!kt || !kt->root || !query || k <= 0)) return;
    if (!dist_fn) dist_fn = spatial_kdtree_dist_sq_euclidean;

    typedef struct { void *data; spatial_num_t dist; } knn_item;

    knn_item *heap = (knn_item*)kt->alloc.malloc(sizeof(knn_item) * (size_t)k, kt->alloc.udata);
    if (SPATIAL_UNLIKELY(!heap)) return;

    int heap_size = 0;
    spatial_num_t worst_dist = SPATIAL_INFINITY;

    #define KNN_SWAP(i, j) do { knn_item t = heap[i]; heap[i] = heap[j]; heap[j] = t; } while(0)
    #define KNN_HEAPIFY_UP(idx) do { \
        int __i = (idx); \
        while (__i > 0 && heap[__i].dist > heap[(__i-1)/2].dist) { \
            KNN_SWAP(__i, (__i-1)/2); __i = (__i-1)/2; \
        } \
    } while(0)
    #define KNN_HEAPIFY_DOWN(idx, sz) do { \
        int __i = (idx); \
        while (2*__i+1 < (sz)) { \
            int j = 2*__i+1; \
            if (j+1 < (sz) && heap[j+1].dist > heap[j].dist) j++; \
            if (heap[__i].dist >= heap[j].dist) break; \
            KNN_SWAP(__i, j); __i = j; \
        } \
    } while(0)

    spatial_kdtree_node *stk[SPATIAL_KDTREE_MAX_DEPTH + 2];
    int sp = 0;
    stk[sp++] = kt->root;

    while (sp > 0) {
        spatial_kdtree_node *node = stk[--sp];

        spatial_num_t node_dist = (spatial_num_t)0.0;
        for (int i = 0; i < SPATIAL_KDTREE_DIMS; i++) {
            if (query[i] < node->min[i]) {
                spatial_num_t d = node->min[i] - query[i];
                node_dist += d * d;
            } else if (query[i] > node->max[i]) {
                spatial_num_t d = query[i] - node->max[i];
                node_dist += d * d;
            }
        }
        if (node_dist >= worst_dist) continue;

        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                spatial_kdtree_item *item = &node->items[i];
                spatial_num_t d = dist_fn(query, item->point, SPATIAL_KDTREE_DIMS, dist_udata);

                if (d < worst_dist) {
                    if (heap_size < k) {
                        heap[heap_size].data = item->data;
                        heap[heap_size].dist = d;
                        KNN_HEAPIFY_UP(heap_size);
                        heap_size++;
                        if (heap_size == k) worst_dist = heap[0].dist;
                    } else {
                        heap[0].data = item->data;
                        heap[0].dist = d;
                        KNN_HEAPIFY_DOWN(0, k);
                        worst_dist = heap[0].dist;
                    }
                }
            }
        } else {
            spatial_num_t diff = query[node->split_axis] - node->split_pos;
            if (diff <= 0) {
                if (node->right) stk[sp++] = node->right;
                if (node->left)  stk[sp++] = node->left;
            } else {
                if (node->left)  stk[sp++] = node->left;
                if (node->right) stk[sp++] = node->right;
            }
        }
    }

    for (int i = 0; i < heap_size && i < k; i++) {
        results[i] = heap[i].data;
        if (distances) distances[i] = heap[i].dist;
    }

    kt->alloc.free(heap, kt->alloc.udata);

    #undef KNN_SWAP
    #undef KNN_HEAPIFY_UP
    #undef KNN_HEAPIFY_DOWN
}

/* Delete (simple for now — removes from leaf only) */
SPATIAL_INLINE bool spatial_kdtree_delete(spatial_kdtree *kt,
                                           const spatial_num_t *min,
                                           const spatial_num_t *max,
                                           spatial_data_t data)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!kt || !kt->root)) return false;

    if (kt->root->is_leaf) {
        for (int i = 0; i < kt->root->item_count; i++) {
            if (kt->root->items[i].data == data) {
                if (kt->callbacks.free) {
                    kt->callbacks.free(kt->root->items[i].data, kt->udata);
                }
                kt->root->items[i] = kt->root->items[kt->root->item_count - 1];
                kt->root->item_count--;
                kt->count--;
                return true;
            }
        }
    }
    return false;
}

typedef bool (*spatial_kdtree_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_kdtree_delete_with_comparator(spatial_kdtree *kt,
                                                           const spatial_num_t *min,
                                                           const spatial_num_t *max,
                                                           spatial_kdtree_cmp_fn cmp,
                                                           void *cmp_udata)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!kt || !kt->root)) return false;

    if (kt->root->is_leaf) {
        for (int i = 0; i < kt->root->item_count; i++) {
            if (cmp(kt->root->items[i].data, cmp_udata)) {
                if (kt->callbacks.free) {
                    kt->callbacks.free(kt->root->items[i].data, kt->udata);
                }
                kt->root->items[i] = kt->root->items[kt->root->item_count - 1];
                kt->root->item_count--;
                kt->count--;
                return true;
            }
        }
    }
    return false;
}

SPATIAL_INLINE void spatial_kdtree_free(spatial_kdtree *kt) {
    if (SPATIAL_UNLIKELY(!kt)) return;
    if (spatial_refcount_decrement(&kt->refcount) > 0) return;
    spatial_kdtree_node_free(kt->root, &kt->alloc, &kt->callbacks, kt->udata);
    kt->alloc.free(kt, kt->alloc.udata);
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_KDTREE_H */

