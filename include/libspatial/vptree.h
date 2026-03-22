/**
 * @file   vptree.h
 * @brief  Vantage-Point Tree for metric space indexing
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 */

#ifndef SPATIAL_VPTREE_H
#define SPATIAL_VPTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"
#include <math.h>

#ifndef SPATIAL_VPTREE_DIMS
    #define SPATIAL_VPTREE_DIMS 3
#endif

#ifndef SPATIAL_VPTREE_MAX_DEPTH
    #define SPATIAL_VPTREE_MAX_DEPTH 32
#endif

#ifndef SPATIAL_VPTREE_LEAF_SIZE
    #define SPATIAL_VPTREE_LEAF_SIZE 8
#endif

_Static_assert(SPATIAL_VPTREE_DIMS > 0, "VP-Tree requires positive dimensionality");
_Static_assert(SPATIAL_VPTREE_LEAF_SIZE >= 1, "Leaf size must be at least 1");

/* Distance function type */
typedef spatial_num_t (*spatial_vptree_dist_fn)(const spatial_num_t *a,
                                                 const spatial_num_t *b,
                                                 int dims,
                                                 void *udata);

/* Forward declarations */
struct spatial_vptree_node;

/* Item structure */
typedef struct spatial_vptree_item {
    spatial_num_t point[SPATIAL_VPTREE_DIMS];
    spatial_data_t data;
} spatial_vptree_item;

/* Node structure - using union for space efficiency */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_vptree_node {
    /* Bounds for pruning */
    spatial_num_t min[SPATIAL_VPTREE_DIMS];
    spatial_num_t max[SPATIAL_VPTREE_DIMS];

    /* Vantage point info (valid when is_leaf=false) */
    spatial_num_t vantage_point[SPATIAL_VPTREE_DIMS];
    spatial_data_t vantage_data;
    spatial_num_t mu;

    /* Union: either children (internal) or items (leaf) */
    union {
        struct {
            struct spatial_vptree_node *SPATIAL_RESTRICT left;
            struct spatial_vptree_node *SPATIAL_RESTRICT right;
        };
        struct {
            spatial_vptree_item *SPATIAL_RESTRICT items;
            int item_count;
        };
    };

    int depth;
    bool is_leaf;
} spatial_vptree_node;

/* Main VP-Tree structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_vptree {
    spatial_vptree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    spatial_vptree_dist_fn dist_fn;
    void *dist_udata;
    int count;
    int max_depth;
    int leaf_size;
    bool relaxed_atomics;
} spatial_vptree;

/* Iterator callback */
typedef bool (*spatial_vptree_iter_fn)(const spatial_num_t *min,
                                        const spatial_num_t *max,
                                        spatial_data_t data,
                                        void *udata);

/* Default distance function (Euclidean squared) */
SPATIAL_INLINE SPATIAL_CONST spatial_num_t spatial_vptree_dist_sq_euclidean(
    const spatial_num_t *a,
    const spatial_num_t *b,
    int dims,
    void *udata)
{
    (void)udata;
    spatial_num_t sum = (spatial_num_t)0.0;
    for (int i = 0; i < dims; i++) {
        spatial_num_t diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sum;
}

/* ============================================================================
 * Internal Helpers
 * ============================================================================ */

SPATIAL_INLINE void spatial_vptree_bounds_init(spatial_vptree_node *node) {
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        node->min[i] = SPATIAL_INFINITY;
        node->max[i] = -SPATIAL_INFINITY;
    }
}

SPATIAL_INLINE void spatial_vptree_bounds_include(spatial_vptree_node *node, 
                                                   const spatial_num_t *point) {
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        if (point[i] < node->min[i]) node->min[i] = point[i];
        if (point[i] > node->max[i]) node->max[i] = point[i];
    }
}

SPATIAL_INLINE spatial_vptree_node* spatial_vptree_node_new_leaf(int depth, 
                                                                  const spatial_allocator *alloc) {
    spatial_vptree_node *node = (spatial_vptree_node*)alloc->malloc(
        sizeof(spatial_vptree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    memset(node, 0, sizeof(spatial_vptree_node));
    node->is_leaf = true;
    node->depth = depth;
    spatial_vptree_bounds_init(node);
    return node;
}

SPATIAL_INLINE void spatial_vptree_node_free(spatial_vptree_node *node,
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
        if (cb && cb->free) {
            cb->free(node->vantage_data, udata);
        }
        spatial_vptree_node_free(node->left, alloc, cb, udata);
        spatial_vptree_node_free(node->right, alloc, cb, udata);
    }
    alloc->free(node, alloc->udata);
}

/* Quickselect for finding median */
SPATIAL_INLINE spatial_num_t spatial_vptree_quickselect(spatial_num_t *arr,
                                                         int left,
                                                         int right,
                                                         int k)
{
    while (left < right) {
        spatial_num_t pivot = arr[(left + right) / 2];
        int i = left, j = right;
        
        while (i <= j) {
            while (arr[i] < pivot) i++;
            while (arr[j] > pivot) j--;
            if (i <= j) {
                spatial_num_t tmp = arr[i];
                arr[i] = arr[j];
                arr[j] = tmp;
                i++;
                j--;
            }
        }
        
        if (k <= j) right = j;
        else if (k >= i) left = i;
        else break;
    }
    
    return arr[k];
}

/* Choose vantage point - furthest from centroid */
SPATIAL_INLINE int spatial_vptree_choose_vantage(spatial_vptree_item *items, int count) {
    if (count <= 1) return 0;

    spatial_num_t centroid[SPATIAL_VPTREE_DIMS];
    for (int d = 0; d < SPATIAL_VPTREE_DIMS; d++) {
        centroid[d] = (spatial_num_t)0.0;
        for (int i = 0; i < count; i++) {
            centroid[d] += items[i].point[d];
        }
        centroid[d] /= (spatial_num_t)count;
    }

    int best_idx = 0;
    spatial_num_t best_dist = (spatial_num_t)0.0;

    for (int i = 0; i < count; i++) {
        spatial_num_t dist = (spatial_num_t)0.0;
        for (int d = 0; d < SPATIAL_VPTREE_DIMS; d++) {
            spatial_num_t diff = items[i].point[d] - centroid[d];
            dist += diff * diff;
        }
        if (dist > best_dist) {
            best_dist = dist;
            best_idx = i;
        }
    }

    return best_idx;
}

SPATIAL_INLINE bool spatial_vptree_node_split(spatial_vptree_node *node,
                                               int max_depth,
                                               const spatial_allocator *alloc,
                                               spatial_vptree_dist_fn dist_fn,
                                               void *dist_udata)
{
    if (SPATIAL_UNLIKELY(!node || !node->is_leaf)) return false;
    if (node->depth >= max_depth) return false;
    if (node->item_count <= 1) return false;

    int count = node->item_count;

    /* Choose vantage point (swap to front) */
    int vp_idx = spatial_vptree_choose_vantage(node->items, count);
    if (vp_idx != 0) {
        spatial_vptree_item tmp = node->items[0];
        node->items[0] = node->items[vp_idx];
        node->items[vp_idx] = tmp;
    }

    /* Compute all distances from vantage point */
    spatial_num_t *distances = (spatial_num_t*)alloc->malloc(
        sizeof(spatial_num_t) * (size_t)(count - 1), alloc->udata);
    if (SPATIAL_UNLIKELY(!distances)) return false;

    for (int i = 1; i < count; i++) {
        distances[i - 1] = dist_fn(node->items[0].point, node->items[i].point, 
                                   SPATIAL_VPTREE_DIMS, dist_udata);
    }

    /* Find median distance */
    spatial_num_t mu;
    if (count > 1) {
        mu = spatial_vptree_quickselect(distances, 0, count - 2, (count - 1) / 2);
    } else {
        mu = (spatial_num_t)0.0;
    }

    /* Count inside/outside */
    int inside_count = 0;
    for (int i = 0; i < count - 1; i++) {
        if (distances[i] <= mu) inside_count++;
    }

    /* Handle degenerate case - ensure both sides get at least one item if possible */
    if (inside_count == 0 && count > 1) {
        inside_count = 1;
        /* Find closest point for inside */
        spatial_num_t min_dist = distances[0];
        int min_idx = 0;
        for (int i = 1; i < count - 1; i++) {
            if (distances[i] < min_dist) {
                min_dist = distances[i];
                min_idx = i;
            }
        }
        mu = distances[min_idx];
    } else if (inside_count >= count - 1 && count > 1) {
        inside_count = count - 2;
        if (inside_count < 1) inside_count = 1;
        /* Recompute mu */
        spatial_num_t *sorted = (spatial_num_t*)alloc->malloc(
            sizeof(spatial_num_t) * (size_t)(count - 1), alloc->udata);
        if (sorted) {
            memcpy(sorted, distances, sizeof(spatial_num_t) * (size_t)(count - 1));
            for (int i = 1; i < count - 1; i++) {
                spatial_num_t key = sorted[i];
                int j = i - 1;
                while (j >= 0 && sorted[j] > key) {
                    sorted[j + 1] = sorted[j];
                    j--;
                }
                sorted[j + 1] = key;
            }
            if (inside_count > 0 && inside_count <= count - 1) {
                mu = sorted[inside_count - 1];
            }
            alloc->free(sorted, alloc->udata);
        }
    }

    int outside_count = count - 1 - inside_count;

    /* Allocate child item arrays */
    spatial_vptree_item *inside_items = NULL;
    spatial_vptree_item *outside_items = NULL;
    
    if (inside_count > 0) {
        inside_items = (spatial_vptree_item*)alloc->malloc(
            sizeof(spatial_vptree_item) * (size_t)inside_count, alloc->udata);
    }
    if (outside_count > 0) {
        outside_items = (spatial_vptree_item*)alloc->malloc(
            sizeof(spatial_vptree_item) * (size_t)outside_count, alloc->udata);
    }

    if ((inside_count > 0 && !inside_items) || (outside_count > 0 && !outside_items)) {
        alloc->free(inside_items, alloc->udata);
        alloc->free(outside_items, alloc->udata);
        alloc->free(distances, alloc->udata);
        return false;
    }

    /* Partition items */
    int ii = 0, oi = 0;
    for (int i = 1; i < count; i++) {
        spatial_num_t d = distances[i - 1];
        if (d <= mu && ii < inside_count) {
            inside_items[ii++] = node->items[i];
        } else {
            outside_items[oi++] = node->items[i];
        }
    }
    alloc->free(distances, alloc->udata);

    /* Create child nodes */
    spatial_vptree_node *inside_node = spatial_vptree_node_new_leaf(node->depth + 1, alloc);
    spatial_vptree_node *outside_node = spatial_vptree_node_new_leaf(node->depth + 1, alloc);

    if (SPATIAL_UNLIKELY(!inside_node || !outside_node)) {
        alloc->free(inside_items, alloc->udata);
        alloc->free(outside_items, alloc->udata);
        spatial_vptree_node_free(inside_node, alloc, NULL, NULL);
        spatial_vptree_node_free(outside_node, alloc, NULL, NULL);
        return false;
    }

    inside_node->items = inside_items;
    inside_node->item_count = inside_count;
    for (int i = 0; i < inside_count; i++) {
        spatial_vptree_bounds_include(inside_node, inside_items[i].point);
    }

    outside_node->items = outside_items;
    outside_node->item_count = outside_count;
    for (int i = 0; i < outside_count; i++) {
        spatial_vptree_bounds_include(outside_node, outside_items[i].point);
    }

    spatial_vptree_item *old_items_array = node->items;

    /* Capture vantage point data before we lose access via the union */
    memcpy(node->vantage_point, old_items_array[0].point, 
           sizeof(spatial_num_t) * SPATIAL_VPTREE_DIMS);
    node->vantage_data = old_items_array[0].data;

    /* Now convert to internal node - this overwrites node->items via union! */
    node->mu = mu;
    node->is_leaf = false;
    node->left = inside_node;    /* Overwrites node->items */
    node->right = outside_node;  /* Overwrites node->item_count */

    /* Update bounds from children */
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        node->min[i] = (inside_node->min[i] < outside_node->min[i]) ? 
                       inside_node->min[i] : outside_node->min[i];
        node->max[i] = (inside_node->max[i] > outside_node->max[i]) ? 
                       inside_node->max[i] : outside_node->max[i];
    }

    /* NOW free the old items array using the cached pointer */
    alloc->free(old_items_array, alloc->udata);
    
    return true;
}

/* ============================================================================
 * Public API
 * ============================================================================ */

SPATIAL_INLINE spatial_vptree* spatial_vptree_new_with_allocator(const spatial_allocator *alloc);

SPATIAL_INLINE spatial_vptree* spatial_vptree_new(void) {
    return spatial_vptree_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_vptree* spatial_vptree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();
    
    spatial_vptree *vp = (spatial_vptree*)alloc->malloc(sizeof(spatial_vptree), alloc->udata);
    if (SPATIAL_UNLIKELY(!vp)) return NULL;
    
    memset(vp, 0, sizeof(spatial_vptree));
    spatial_refcount_init(&vp->refcount);
    vp->alloc = *alloc;
    vp->dist_fn = spatial_vptree_dist_sq_euclidean;
    vp->max_depth = SPATIAL_VPTREE_MAX_DEPTH;
    vp->leaf_size = SPATIAL_VPTREE_LEAF_SIZE;
    vp->relaxed_atomics = true;
    vp->root = NULL;
    
    return vp;
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_vptree* spatial_vptree_new_with_memento_arena(memento_arena_t *arena) {
    spatial_allocator alloc = spatial_memento_arena_allocator(arena);
    return spatial_vptree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_vptree* spatial_vptree_new_with_memento_pool(memento_pool_t *pool) {
    spatial_allocator alloc = spatial_memento_pool_allocator(pool);
    return spatial_vptree_new_with_allocator(&alloc);
}
#endif

SPATIAL_INLINE spatial_vptree* spatial_vptree_new_with_arena(spatial_arena_t *arena) {
    spatial_allocator alloc = spatial_arena_allocator(arena);
    return spatial_vptree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_vptree* spatial_vptree_new_with_pool(spatial_pool_t *pool) {
    spatial_allocator alloc = spatial_pool_allocator(pool);
    return spatial_vptree_new_with_allocator(&alloc);
}

SPATIAL_INLINE void spatial_vptree_free(spatial_vptree *vp) {
    if (SPATIAL_UNLIKELY(!vp)) return;
    
    if (spatial_refcount_decrement(&vp->refcount) > 0) return;
    
    spatial_vptree_node_free(vp->root, &vp->alloc, &vp->callbacks, vp->udata);
    vp->alloc.free(vp, vp->alloc.udata);
}

SPATIAL_INLINE spatial_vptree* spatial_vptree_clone(const spatial_vptree *vp) {
    if (SPATIAL_UNLIKELY(!vp)) return NULL;
    
    spatial_refcount_increment(&vp->refcount);
    union { const spatial_vptree *in; spatial_vptree *out; } pun; 
    pun.in = vp; 
    return pun.out;
}

SPATIAL_INLINE void spatial_vptree_set_item_callbacks(spatial_vptree *vp,
                                                       const spatial_item_callbacks *cb)
{
    if (vp) vp->callbacks = *cb;
}

SPATIAL_INLINE void spatial_vptree_set_udata(spatial_vptree *vp, void *udata) {
    if (vp) vp->udata = udata;
}

SPATIAL_INLINE void spatial_vptree_set_distance(spatial_vptree *vp,
                                                 spatial_vptree_dist_fn fn,
                                                 void *dist_udata)
{
    if (vp) {
        vp->dist_fn = fn ? fn : spatial_vptree_dist_sq_euclidean;
        vp->dist_udata = dist_udata;
    }
}

SPATIAL_INLINE void spatial_vptree_opt_relaxed_atomics(spatial_vptree *vp, bool enable) {
    if (vp) vp->relaxed_atomics = enable;
}

/* ============================================================================
 * INCREMENTAL INSERT - with safety checks
 * ============================================================================ */

SPATIAL_INLINE bool spatial_vptree_insert(spatial_vptree *vp,
                                           const spatial_num_t *min,
                                           const spatial_num_t *max,
                                           spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!vp)) return false;
    
    /* Compute point from min/max */
    spatial_num_t point[SPATIAL_VPTREE_DIMS];
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        point[i] = (min[i] + max[i]) * (spatial_num_t)0.5;
    }
    
    /* Create root if needed */
    if (SPATIAL_UNLIKELY(!vp->root)) {
        vp->root = spatial_vptree_node_new_leaf(0, &vp->alloc);
        if (SPATIAL_UNLIKELY(!vp->root)) return false;
        for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
            vp->root->min[i] = point[i];
            vp->root->max[i] = point[i];
        }
    }
    
    /* Navigate to leaf with safety checks */
    spatial_vptree_node *node = vp->root;
    int depth = 0;
    
    while (node && !node->is_leaf) {
        /* Safety: prevent infinite loop */
        if (SPATIAL_UNLIKELY(depth++ > vp->max_depth * 2)) {
            return false;
        }
        
        /* Update bounds */
        spatial_vptree_bounds_include(node, point);
        
        spatial_num_t d = vp->dist_fn(point, node->vantage_point, 
                                      SPATIAL_VPTREE_DIMS, vp->dist_udata);
        
        /* Check for NULL children before navigating */
        if (d <= node->mu) {
            if (SPATIAL_UNLIKELY(!node->left)) {
                node->left = spatial_vptree_node_new_leaf(node->depth + 1, &vp->alloc);
                if (!node->left) return false;
            }
            node = node->left;
        } else {
            if (SPATIAL_UNLIKELY(!node->right)) {
                node->right = spatial_vptree_node_new_leaf(node->depth + 1, &vp->alloc);
                if (!node->right) return false;
            }
            node = node->right;
        }
    }
    
    /* Validate we reached a leaf */
    if (SPATIAL_UNLIKELY(!node || !node->is_leaf)) {
        return false;
    }
    
    /* Update leaf bounds */
    spatial_vptree_bounds_include(node, point);
    
    /* Grow item array */
    int new_count = node->item_count + 1;
    spatial_vptree_item *new_items = (spatial_vptree_item*)vp->alloc.malloc(
        sizeof(spatial_vptree_item) * (size_t)new_count, vp->alloc.udata);
    if (SPATIAL_UNLIKELY(!new_items)) return false;
    
    if (node->items) {
        memcpy(new_items, node->items, 
               sizeof(spatial_vptree_item) * (size_t)node->item_count);
        vp->alloc.free(node->items, vp->alloc.udata);
    }
    
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        new_items[node->item_count].point[i] = point[i];
    }
    new_items[node->item_count].data = data;
    
    node->items = new_items;
    node->item_count = new_count;
    vp->count++;
    
    /* Split if overflow */
    if (node->item_count > vp->leaf_size && node->depth < vp->max_depth) {
        (void)spatial_vptree_node_split(node, vp->max_depth, &vp->alloc, 
                                        vp->dist_fn, vp->dist_udata);
    }
    
    return true;
}

/* ============================================================================
 * Bulk build for initial data loading
 * ============================================================================ */

SPATIAL_INLINE spatial_vptree_node* spatial_vptree_build_recursive(
    spatial_vptree_item *items,
    int count,
    int depth,
    int max_depth,
    int leaf_size,
    spatial_vptree_dist_fn dist_fn,
    void *dist_udata,
    const spatial_allocator *alloc)
{
    if (count <= leaf_size || depth >= max_depth) {
        spatial_vptree_node *leaf = spatial_vptree_node_new_leaf(depth, alloc);
        if (SPATIAL_UNLIKELY(!leaf)) return NULL;
        
        leaf->items = (spatial_vptree_item*)alloc->malloc(
            sizeof(spatial_vptree_item) * (size_t)count, alloc->udata);
        if (SPATIAL_UNLIKELY(!leaf->items)) {
            alloc->free(leaf, alloc->udata);
            return NULL;
        }
        
        memcpy(leaf->items, items, sizeof(spatial_vptree_item) * (size_t)count);
        leaf->item_count = count;
        
        for (int i = 0; i < count; i++) {
            spatial_vptree_bounds_include(leaf, items[i].point);
        }
        
        return leaf;
    }
    
    /* Choose vantage point */
    int vp_idx = spatial_vptree_choose_vantage(items, count);
    if (vp_idx != 0) {
        spatial_vptree_item tmp = items[0];
        items[0] = items[vp_idx];
        items[vp_idx] = tmp;
    }
    
    spatial_vptree_node *node = (spatial_vptree_node*)alloc->malloc(
        sizeof(spatial_vptree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    memset(node, 0, sizeof(spatial_vptree_node));
    node->is_leaf = false;
    node->depth = depth;
    spatial_vptree_bounds_init(node);
    
    memcpy(node->vantage_point, items[0].point, 
           sizeof(spatial_num_t) * SPATIAL_VPTREE_DIMS);
    node->vantage_data = items[0].data;
    
    /* Compute distances */
    spatial_num_t *distances = (spatial_num_t*)alloc->malloc(
        sizeof(spatial_num_t) * (size_t)(count - 1), alloc->udata);
    if (SPATIAL_UNLIKELY(!distances)) {
        alloc->free(node, alloc->udata);
        return NULL;
    }
    
    for (int i = 1; i < count; i++) {
        distances[i - 1] = dist_fn(node->vantage_point, items[i].point, 
                                   SPATIAL_VPTREE_DIMS, dist_udata);
    }
    
    /* Find median */
    spatial_num_t mu;
    if (count > 1) {
        mu = spatial_vptree_quickselect(distances, 0, count - 2, (count - 1) / 2);
    } else {
        mu = (spatial_num_t)0.0;
    }
    node->mu = mu;
    
    /* Count and partition */
    int inside_count = 0;
    for (int i = 0; i < count - 1; i++) {
        if (distances[i] <= mu) inside_count++;
    }
    
    /* Handle degenerate case */
    if (inside_count == 0 || inside_count >= count - 1) {
        inside_count = (count - 1) / 2;
        if (inside_count < 1 && count > 1) inside_count = 1;
        
        spatial_num_t *sorted = (spatial_num_t*)alloc->malloc(
            sizeof(spatial_num_t) * (size_t)(count - 1), alloc->udata);
        if (sorted) {
            memcpy(sorted, distances, sizeof(spatial_num_t) * (size_t)(count - 1));
            for (int i = 1; i < count - 1; i++) {
                spatial_num_t key = sorted[i];
                int j = i - 1;
                while (j >= 0 && sorted[j] > key) {
                    sorted[j + 1] = sorted[j];
                    j--;
                }
                sorted[j + 1] = key;
            }
            if (inside_count > 0 && inside_count <= count - 1) {
                mu = sorted[inside_count - 1];
            }
            alloc->free(sorted, alloc->udata);
        }
        node->mu = mu;
    }
    
    int outside_count = count - 1 - inside_count;
    
    spatial_vptree_item *inside_items = NULL;
    spatial_vptree_item *outside_items = NULL;
    
    if (inside_count > 0) {
        inside_items = (spatial_vptree_item*)alloc->malloc(
            sizeof(spatial_vptree_item) * (size_t)inside_count, alloc->udata);
    }
    if (outside_count > 0) {
        outside_items = (spatial_vptree_item*)alloc->malloc(
            sizeof(spatial_vptree_item) * (size_t)outside_count, alloc->udata);
    }
    
    if ((inside_count > 0 && !inside_items) || (outside_count > 0 && !outside_items)) {
        alloc->free(inside_items, alloc->udata);
        alloc->free(outside_items, alloc->udata);
        alloc->free(distances, alloc->udata);
        alloc->free(node, alloc->udata);
        return NULL;
    }
    
    int ii = 0, oi = 0;
    for (int i = 1; i < count; i++) {
        spatial_num_t d = distances[i - 1];
        if (d <= mu && ii < inside_count) {
            inside_items[ii++] = items[i];
        } else {
            outside_items[oi++] = items[i];
        }
    }
    
    alloc->free(distances, alloc->udata);
    
    node->left = spatial_vptree_build_recursive(inside_items, inside_count, depth + 1,
                                                   max_depth, leaf_size, dist_fn, dist_udata, alloc);
    node->right = spatial_vptree_build_recursive(outside_items, outside_count, depth + 1,
                                                    max_depth, leaf_size, dist_fn, dist_udata, alloc);
    
    alloc->free(inside_items, alloc->udata);
    alloc->free(outside_items, alloc->udata);
    
    if (SPATIAL_UNLIKELY(!node->left || !node->right)) {
        spatial_vptree_node_free(node->left, alloc, NULL, NULL);
        spatial_vptree_node_free(node->right, alloc, NULL, NULL);
        alloc->free(node, alloc->udata);
        return NULL;
    }
    
    /* Compute bounds from children */
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        node->min[i] = (node->left->min[i] < node->right->min[i]) ? 
                       node->left->min[i] : node->right->min[i];
        node->max[i] = (node->left->max[i] > node->right->max[i]) ? 
                       node->left->max[i] : node->right->max[i];
    }
    
    return node;
}

/* ============================================================================
 * Search operations
 * ============================================================================ */

SPATIAL_INLINE void spatial_vptree_search(const spatial_vptree *vp,
                                           const spatial_num_t *qmin,
                                           const spatial_num_t *qmax,
                                           spatial_vptree_iter_fn iter,
                                           void *udata)
{
    if (SPATIAL_UNLIKELY(!vp || !vp->root || !iter)) return;
    
    spatial_num_t query[SPATIAL_VPTREE_DIMS];
    spatial_num_t radius = (spatial_num_t)0.0;
    
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        query[i] = (qmin[i] + qmax[i]) * (spatial_num_t)0.5;
        spatial_num_t half_size = (qmax[i] - qmin[i]) * (spatial_num_t)0.5;
        radius += half_size * half_size;
    }
    radius = (spatial_num_t)sqrt((double)radius);
    spatial_num_t radius_sq = radius * radius;
    
    spatial_vptree_node *stk[SPATIAL_VPTREE_MAX_DEPTH + 2];
    int sp = 0;
    stk[sp++] = vp->root;

    while (sp > 0) {
        spatial_vptree_node *node = stk[--sp];
        if (!node) continue;

        /* Bounds pruning */
        spatial_num_t node_dist_sq = (spatial_num_t)0.0;
        for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
            if (query[i] < node->min[i]) {
                spatial_num_t d = node->min[i] - query[i];
                node_dist_sq += d * d;
            } else if (query[i] > node->max[i]) {
                spatial_num_t d = query[i] - node->max[i];
                node_dist_sq += d * d;
            }
        }
        if (node_dist_sq > radius_sq) continue;

        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                spatial_num_t d = vp->dist_fn(query, node->items[i].point,
                                              SPATIAL_VPTREE_DIMS, vp->dist_udata);
                if (d <= radius_sq) {
                    if (!iter(node->items[i].point, node->items[i].point,
                              node->items[i].data, udata)) return;
                }
            }
        } else {
            spatial_num_t d = vp->dist_fn(query, node->vantage_point,
                                          SPATIAL_VPTREE_DIMS, vp->dist_udata);

            if (d <= radius_sq) {
                if (!iter(node->vantage_point, node->vantage_point,
                          node->vantage_data, udata)) return;
            }

            spatial_num_t diff = d - node->mu;
            if (diff <= radius  && node->left)  stk[sp++] = node->left;
            if (diff >= -radius && node->right) stk[sp++] = node->right;
        }
    }
}

/* Scan all items */
SPATIAL_INLINE void spatial_vptree_scan(const spatial_vptree *vp,
                                         spatial_vptree_iter_fn iter,
                                         void *udata)
{
    if (SPATIAL_UNLIKELY(!vp || !vp->root || !iter)) return;

    spatial_vptree_node *stk[SPATIAL_VPTREE_MAX_DEPTH + 2];
    int sp = 0;
    stk[sp++] = vp->root;

    while (sp > 0) {
        spatial_vptree_node *node = stk[--sp];
        if (!node) continue;

        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                if (!iter(node->items[i].point, node->items[i].point,
                          node->items[i].data, udata)) return;
            }
        } else {
            if (!iter(node->vantage_point, node->vantage_point,
                      node->vantage_data, udata)) return;
            if (node->left)  stk[sp++] = node->left;
            if (node->right) stk[sp++] = node->right;
        }
    }
}

SPATIAL_INLINE int spatial_vptree_count(const spatial_vptree *vp) {
    return vp ? vp->count : 0;
}

/* ============================================================================
 * k-NN search
 * ============================================================================ */

SPATIAL_INLINE void spatial_vptree_nearest(
    const spatial_vptree *vp,
    const spatial_num_t *query,
    int k,
    void **results,
    spatial_num_t *distances)
{
    if (SPATIAL_UNLIKELY(!vp || !vp->root || !query || k <= 0)) return;
    
    typedef struct {
        void *data;
        spatial_num_t dist;
    } knn_item;
    
    knn_item *heap = (knn_item*)vp->alloc.malloc(sizeof(knn_item) * (size_t)k, vp->alloc.udata);
    if (SPATIAL_UNLIKELY(!heap)) return;
    
    int heap_size = 0;
    spatial_num_t worst_dist = SPATIAL_INFINITY;
    
    #define VP_KNN_SWAP(i, j) do { knn_item t = heap[i]; heap[i] = heap[j]; heap[j] = t; } while(0)
    
    #define VP_KNN_HEAPIFY_UP(idx) do { \
        int __vp_i = (idx); \
        while (__vp_i > 0 && heap[__vp_i].dist > heap[(__vp_i-1)/2].dist) { \
            VP_KNN_SWAP(__vp_i, (__vp_i-1)/2); \
            __vp_i = (__vp_i-1)/2; \
        } \
    } while(0)
    
    #define VP_KNN_HEAPIFY_DOWN(idx, sz) do { \
        int __vp_i = (idx); \
        while (2*__vp_i+1 < (sz)) { \
            int __vp_j = 2*__vp_i+1; \
            if (__vp_j+1 < (sz) && heap[__vp_j+1].dist > heap[__vp_j].dist) __vp_j++; \
            if (heap[__vp_i].dist >= heap[__vp_j].dist) break; \
            VP_KNN_SWAP(__vp_i, __vp_j); \
            __vp_i = __vp_j; \
        } \
    } while(0)
    
    spatial_vptree_node *stk[SPATIAL_VPTREE_MAX_DEPTH + 2];
    int sp = 0;

    /* Push root with bounds check */
    spatial_num_t root_dist_sq = (spatial_num_t)0.0;
    for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
        if (query[i] < vp->root->min[i]) {
            spatial_num_t d = vp->root->min[i] - query[i];
            root_dist_sq += d * d;
        } else if (query[i] > vp->root->max[i]) {
            spatial_num_t d = query[i] - vp->root->max[i];
            root_dist_sq += d * d;
        }
    }

    if (root_dist_sq < worst_dist) {
        stk[sp++] = vp->root;
    }

    while (sp > 0) {
        spatial_vptree_node *node = stk[--sp];
        if (!node) continue;
        
        /* Bounds pruning */
        spatial_num_t node_dist_sq = (spatial_num_t)0.0;
        for (int i = 0; i < SPATIAL_VPTREE_DIMS; i++) {
            if (query[i] < node->min[i]) {
                spatial_num_t d = node->min[i] - query[i];
                node_dist_sq += d * d;
            } else if (query[i] > node->max[i]) {
                spatial_num_t d = query[i] - node->max[i];
                node_dist_sq += d * d;
            }
        }
        if (node_dist_sq >= worst_dist) continue;
        
        if (node->is_leaf) {
            for (int i = 0; i < node->item_count; i++) {
                spatial_num_t d = vp->dist_fn(query, node->items[i].point, 
                                              SPATIAL_VPTREE_DIMS, vp->dist_udata);
                
                if (d < worst_dist) {
                    if (heap_size < k) {
                        heap[heap_size].data = node->items[i].data;
                        heap[heap_size].dist = d;
                        VP_KNN_HEAPIFY_UP(heap_size);
                        heap_size++;
                        if (heap_size == k) worst_dist = heap[0].dist;
                    } else {
                        heap[0].data = node->items[i].data;
                        heap[0].dist = d;
                        VP_KNN_HEAPIFY_DOWN(0, k);
                        worst_dist = heap[0].dist;
                    }
                }
            }
        } else {
            spatial_num_t d = vp->dist_fn(query, node->vantage_point, 
                                          SPATIAL_VPTREE_DIMS, vp->dist_udata);
            
            if (d < worst_dist) {
                if (heap_size < k) {
                    heap[heap_size].data = node->vantage_data;
                    heap[heap_size].dist = d;
                    VP_KNN_HEAPIFY_UP(heap_size);
                    heap_size++;
                    if (heap_size == k) worst_dist = heap[0].dist;
                } else {
                    heap[0].data = node->vantage_data;
                    heap[0].dist = d;
                    VP_KNN_HEAPIFY_DOWN(0, k);
                    worst_dist = heap[0].dist;
                }
            }
            
            /* Triangle inequality pruning */
            spatial_num_t diff = d - node->mu;
            spatial_num_t diff_sq = diff * diff;
            
            /* Visit closer side first for better pruning */
            if (diff <= 0) {
                /* Left/inside is closer */
                if (node->right && (diff_sq < worst_dist || heap_size < k))
                    stk[sp++] = node->right;
                if (node->left) stk[sp++] = node->left;
            } else {
                /* Right/outside is closer */
                if (node->left && (diff_sq < worst_dist || heap_size < k))
                    stk[sp++] = node->left;
                if (node->right) stk[sp++] = node->right;
            }
        }
    }
    
    for (int i = 0; i < heap_size && i < k; i++) {
        results[i] = heap[i].data;
        if (distances) distances[i] = heap[i].dist;
    }
    
    vp->alloc.free(heap, vp->alloc.udata);
    
    #undef VP_KNN_SWAP
    #undef VP_KNN_HEAPIFY_UP
    #undef VP_KNN_HEAPIFY_DOWN
}

/* ============================================================================
 * Delete operations
 * ============================================================================ */

SPATIAL_INLINE bool spatial_vptree_delete(spatial_vptree *vp,
                                           const spatial_num_t *min,
                                           const spatial_num_t *max,
                                           spatial_data_t data)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!vp || !vp->root)) return false;
    
    /* For now, only support deletion from root leaf */
    if (vp->root->is_leaf) {
        for (int i = 0; i < vp->root->item_count; i++) {
            if (vp->root->items[i].data == data) {
                if (vp->callbacks.free) {
                    vp->callbacks.free(vp->root->items[i].data, vp->udata);
                }
                vp->root->items[i] = vp->root->items[vp->root->item_count - 1];
                vp->root->item_count--;
                vp->count--;
                return true;
            }
        }
    }
    
    return false;
}

typedef bool (*spatial_vptree_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_vptree_delete_with_comparator(spatial_vptree *vp,
                                                           const spatial_num_t *min,
                                                           const spatial_num_t *max,
                                                           spatial_vptree_cmp_fn cmp,
                                                           void *cmp_udata)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!vp || !vp->root)) return false;
    
    if (vp->root->is_leaf) {
        for (int i = 0; i < vp->root->item_count; i++) {
            if (cmp(vp->root->items[i].data, cmp_udata)) {
                if (vp->callbacks.free) {
                    vp->callbacks.free(vp->root->items[i].data, vp->udata);
                }
                vp->root->items[i] = vp->root->items[vp->root->item_count - 1];
                vp->root->item_count--;
                vp->count--;
                return true;
            }
        }
    }
    
    return false;
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_VPTREE_H */

