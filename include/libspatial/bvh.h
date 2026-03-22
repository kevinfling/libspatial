/**
 * @file   bvh.h
 * @brief  Bounding Volume Hierarchy with SAH build
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * BVH is a tree structure where each node stores a bounding volume that
 * encloses all its children. Uses Surface Area Heuristic (SAH) for optimal
 * build quality.
 */

#ifndef SPATIAL_BVH_H
#define SPATIAL_BVH_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

#ifndef SPATIAL_BVH_DIMS
    #define SPATIAL_BVH_DIMS 3
#endif

#ifndef SPATIAL_BVH_MAX_DEPTH
    #define SPATIAL_BVH_MAX_DEPTH 32
#endif

#ifndef SPATIAL_BVH_LEAF_SIZE
    #define SPATIAL_BVH_LEAF_SIZE 4
#endif

#ifndef SPATIAL_BVH_TRAVERSAL_COST
    #define SPATIAL_BVH_TRAVERSAL_COST 1.0
#endif

#ifndef SPATIAL_BVH_INTERSECTION_COST
    #define SPATIAL_BVH_INTERSECTION_COST 1.0
#endif

SPATIAL_STATIC_ASSERT(bvh_dims_valid,
    SPATIAL_BVH_DIMS > 0 && SPATIAL_BVH_DIMS <= 8);

/* Forward declarations */
struct spatial_bvh_node;

/* Primitive structure (input geometry) */
typedef struct spatial_bvh_prim {
    spatial_num_t min[SPATIAL_BVH_DIMS];
    spatial_num_t max[SPATIAL_BVH_DIMS];
    spatial_num_t center[SPATIAL_BVH_DIMS];
    spatial_data_t data;
    int id;
} spatial_bvh_prim;

/* Node structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_bvh_node {
    spatial_num_t min[SPATIAL_BVH_DIMS];
    spatial_num_t max[SPATIAL_BVH_DIMS];
    
    union {
        struct {
            int left_child;
            int right_child;
        };
        struct {
            int prim_offset;
            int prim_count;
        };
    };
    
    int flags;  /* bit 0: is_leaf */
} spatial_bvh_node;

/* Main BVH structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_bvh {
    spatial_bvh_prim *prims;
    int prim_count;
    int prim_capacity;
    
    spatial_bvh_node *nodes;
    int node_count;
    int node_capacity;
    
    int *indices;  /* Primitive indices for reordering */
    
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    
    int max_depth;
    int leaf_size;
    bool built;
    bool relaxed_atomics;
} spatial_bvh;

/* Iterator callback */
typedef bool (*spatial_bvh_iter_fn)(const spatial_num_t *min,
                                     const spatial_num_t *max,
                                     spatial_data_t data,
                                     void *udata);

/* Forward declarations for public API */
SPATIAL_INLINE spatial_bvh* spatial_bvh_new_with_allocator(const spatial_allocator *alloc);

/* Ray structure for ray intersection */
typedef struct {
    spatial_num_t origin[SPATIAL_BVH_DIMS];
    spatial_num_t dir[SPATIAL_BVH_DIMS];
    spatial_num_t inv_dir[SPATIAL_BVH_DIMS];
    spatial_num_t tmin;
    spatial_num_t tmax;
} spatial_bvh_ray;

/* Ray hit result */
typedef struct {
    spatial_data_t data;
    spatial_num_t t;
    int prim_id;
} spatial_bvh_hit;

/* ============================================================================
 * Ray-Box Intersection (Slab method)
 * ============================================================================ */

SPATIAL_INLINE bool spatial_bvh_ray_intersect_box(const spatial_bvh_ray *SPATIAL_RESTRICT ray,
                                                   const spatial_num_t *SPATIAL_RESTRICT bmin,
                                                   const spatial_num_t *SPATIAL_RESTRICT bmax,
                                                   spatial_num_t tmin,
                                                   spatial_num_t tmax)
{
    for (int i = 0; i < SPATIAL_BVH_DIMS; i++) {
        spatial_num_t t1 = (bmin[i] - ray->origin[i]) * ray->inv_dir[i];
        spatial_num_t t2 = (bmax[i] - ray->origin[i]) * ray->inv_dir[i];
        
        spatial_num_t tn = spatial_min(t1, t2);
        spatial_num_t tf = spatial_max(t1, t2);
        
        tmin = spatial_max(tmin, tn);
        tmax = spatial_min(tmax, tf);
        
        if (tmin > tmax) return false;
    }
    return tmax >= (spatial_num_t)0.0;
}

/* ============================================================================
 * Internal utilities
 * ============================================================================ */

SPATIAL_INLINE void spatial_bvh_node_init(spatial_bvh_node *node)
{
    spatial_bbox_init(node->min, node->max, SPATIAL_BVH_DIMS);
    node->left_child = -1;
    node->right_child = -1;
    node->flags = 0;
}

SPATIAL_INLINE bool spatial_bvh_node_is_leaf(const spatial_bvh_node *node)
{
    return (node->flags & 1) != 0;
}

SPATIAL_INLINE void spatial_bvh_node_set_leaf(spatial_bvh_node *node, bool is_leaf)
{
    if (is_leaf) node->flags |= 1;
    else node->flags &= ~1;
}

SPATIAL_INLINE void spatial_bvh_compute_bounds(spatial_bvh *bvh,
                                                int start,
                                                int end,
                                                spatial_num_t *SPATIAL_RESTRICT out_min,
                                                spatial_num_t *SPATIAL_RESTRICT out_max)
{
    spatial_bbox_init(out_min, out_max, SPATIAL_BVH_DIMS);
    
    for (int i = start; i < end; i++) {
        spatial_bvh_prim *p = &bvh->prims[bvh->indices[i]];
        for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
            out_min[d] = spatial_min(out_min[d], p->min[d]);
            out_max[d] = spatial_max(out_max[d], p->max[d]);
        }
    }
}

SPATIAL_INLINE spatial_num_t spatial_bvh_compute_sa(const spatial_num_t *min,
                                                     const spatial_num_t *max)
{
    return spatial_bbox_area(min, max, SPATIAL_BVH_DIMS);
}

/* Bucket sort for SAH */
#define SPATIAL_BVH_NUM_BUCKETS 12

typedef struct {
    spatial_num_t min[SPATIAL_BVH_DIMS];
    spatial_num_t max[SPATIAL_BVH_DIMS];
    int count;
} spatial_bvh_bucket;

SPATIAL_INLINE spatial_num_t spatial_bvh_eval_sah(spatial_bvh_bucket *buckets,
                                                   int num_buckets,
                                                   spatial_num_t total_sa)
{
    spatial_num_t min_cost = SPATIAL_INFINITY;
    
    for (int i = 1; i < num_buckets; i++) {
        spatial_num_t left_min[SPATIAL_BVH_DIMS];
        spatial_num_t left_max[SPATIAL_BVH_DIMS];
        spatial_num_t right_min[SPATIAL_BVH_DIMS];
        spatial_num_t right_max[SPATIAL_BVH_DIMS];
        
        spatial_bbox_init(left_min, left_max, SPATIAL_BVH_DIMS);
        spatial_bbox_init(right_min, right_max, SPATIAL_BVH_DIMS);
        
        int left_count = 0, right_count = 0;
        
        for (int j = 0; j < i; j++) {
            if (buckets[j].count > 0) {
                for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                    left_min[d] = spatial_min(left_min[d], buckets[j].min[d]);
                    left_max[d] = spatial_max(left_max[d], buckets[j].max[d]);
                }
                left_count += buckets[j].count;
            }
        }
        
        for (int j = i; j < num_buckets; j++) {
            if (buckets[j].count > 0) {
                for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                    right_min[d] = spatial_min(right_min[d], buckets[j].min[d]);
                    right_max[d] = spatial_max(right_max[d], buckets[j].max[d]);
                }
                right_count += buckets[j].count;
            }
        }
        
        if (left_count == 0 || right_count == 0) continue;
        
        spatial_num_t left_sa = spatial_bvh_compute_sa(left_min, left_max);
        spatial_num_t right_sa = spatial_bvh_compute_sa(right_min, right_max);
        
        spatial_num_t cost = spatial_sah_cost(total_sa, left_sa, right_sa,
                                               left_count, right_count,
                                               SPATIAL_BVH_TRAVERSAL_COST,
                                               SPATIAL_BVH_INTERSECTION_COST);
        
        min_cost = spatial_min(min_cost, cost);
    }
    
    return min_cost;
}

/* ============================================================================
 * Public API
 * ============================================================================ */

SPATIAL_INLINE spatial_bvh* spatial_bvh_new(void) {
    return spatial_bvh_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_bvh* spatial_bvh_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();
    
    spatial_bvh *bvh = (spatial_bvh*)alloc->malloc(sizeof(spatial_bvh), alloc->udata);
    if (SPATIAL_UNLIKELY(!bvh)) return NULL;
    
    memset(bvh, 0, sizeof(spatial_bvh));
    spatial_refcount_init(&bvh->refcount);
    bvh->alloc = *alloc;
    bvh->max_depth = SPATIAL_BVH_MAX_DEPTH;
    bvh->leaf_size = SPATIAL_BVH_LEAF_SIZE;
    bvh->relaxed_atomics = true;
    bvh->built = false;
    
    bvh->prim_capacity = 64;
    bvh->prims = (spatial_bvh_prim*)alloc->malloc(
        sizeof(spatial_bvh_prim) * (size_t)bvh->prim_capacity, alloc->udata);
    if (SPATIAL_UNLIKELY(!bvh->prims)) {
        alloc->free(bvh, alloc->udata);
        return NULL;
    }
    
    return bvh;
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_bvh* spatial_bvh_new_with_memento_arena(memento_arena_t *arena) {
    spatial_allocator alloc = spatial_memento_arena_allocator(arena);
    return spatial_bvh_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_bvh* spatial_bvh_new_with_memento_pool(memento_pool_t *pool) {
    spatial_allocator alloc = spatial_memento_pool_allocator(pool);
    return spatial_bvh_new_with_allocator(&alloc);
}
#endif

SPATIAL_INLINE spatial_bvh* spatial_bvh_new_with_arena(spatial_arena_t *arena) {
    spatial_allocator alloc = spatial_arena_allocator(arena);
    return spatial_bvh_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_bvh* spatial_bvh_new_with_pool(spatial_pool_t *pool) {
    spatial_allocator alloc = spatial_pool_allocator(pool);
    return spatial_bvh_new_with_allocator(&alloc);
}

SPATIAL_INLINE void spatial_bvh_free(spatial_bvh *bvh) {
    if (SPATIAL_UNLIKELY(!bvh)) return;
    
    if (spatial_refcount_decrement(&bvh->refcount) > 0) return;
    
    if (bvh->callbacks.free) {
        for (int i = 0; i < bvh->prim_count; i++) {
            bvh->callbacks.free(bvh->prims[i].data, bvh->udata);
        }
    }
    
    bvh->alloc.free(bvh->prims, bvh->alloc.udata);
    bvh->alloc.free(bvh->nodes, bvh->alloc.udata);
    bvh->alloc.free(bvh->indices, bvh->alloc.udata);
    bvh->alloc.free(bvh, bvh->alloc.udata);
}

SPATIAL_INLINE spatial_bvh* spatial_bvh_clone(const spatial_bvh *bvh) {
    if (SPATIAL_UNLIKELY(!bvh)) return NULL;
    
    spatial_refcount_increment(&bvh->refcount);
    union { const spatial_bvh *in; spatial_bvh *out; } pun; pun.in = bvh; return pun.out;
}

SPATIAL_INLINE void spatial_bvh_set_item_callbacks(spatial_bvh *bvh,
                                                    const spatial_item_callbacks *cb)
{
    if (bvh) bvh->callbacks = *cb;
}

SPATIAL_INLINE void spatial_bvh_set_udata(spatial_bvh *bvh, void *udata) {
    if (bvh) bvh->udata = udata;
}

SPATIAL_INLINE void spatial_bvh_opt_relaxed_atomics(spatial_bvh *bvh, bool enable) {
    if (bvh) bvh->relaxed_atomics = enable;
}

SPATIAL_INLINE bool spatial_bvh_insert(spatial_bvh *bvh,
                                        const spatial_num_t *min,
                                        const spatial_num_t *max,
                                        spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!bvh)) return false;
    
    if (bvh->built) return false;  /* Cannot insert after build */
    
    /* Grow array if needed */
    if (bvh->prim_count >= bvh->prim_capacity) {
        int new_cap = bvh->prim_capacity * 2;
        spatial_bvh_prim *new_prims = (spatial_bvh_prim*)bvh->alloc.malloc(
            sizeof(spatial_bvh_prim) * (size_t)new_cap, bvh->alloc.udata);
        if (SPATIAL_UNLIKELY(!new_prims)) return false;
        
        memcpy(new_prims, bvh->prims, sizeof(spatial_bvh_prim) * (size_t)bvh->prim_count);
        bvh->alloc.free(bvh->prims, bvh->alloc.udata);
        bvh->prims = new_prims;
        bvh->prim_capacity = new_cap;
    }
    
    spatial_bvh_prim *p = &bvh->prims[bvh->prim_count];
    spatial_bbox_copy(min, max, p->min, p->max, SPATIAL_BVH_DIMS);
    
    for (int i = 0; i < SPATIAL_BVH_DIMS; i++) {
        p->center[i] = (min[i] + max[i]) * (spatial_num_t)0.5;
    }
    
    p->data = data;
    p->id = bvh->prim_count;
    
    bvh->prim_count++;
    return true;
}

SPATIAL_INLINE void spatial_bvh_median_split(spatial_bvh *bvh,
                                              int *indices,
                                              int start,
                                              int end,
                                              int axis,
                                              double *scratch_d,
                                              int *scratch_i)
{
    int count = end - start;
    if (count <= 1) return;

    for (int i = 0; i < count; i++) {
        scratch_d[i] = (double)bvh->prims[indices[start + i]].center[axis];
        scratch_i[i] = indices[start + i];
    }

    int mid = count / 2;
    spatial_quickselect_double(scratch_d, scratch_i, 0, count - 1, mid);

    for (int i = 0; i < count; i++) {
        indices[start + i] = scratch_i[i];
    }
}

/* Build BVH using SAH with quickselect median split */
SPATIAL_INLINE int spatial_bvh_build_recursive(spatial_bvh *bvh,
                                                int start,
                                                int end,
                                                int depth,
                                                double *scratch_d,
                                                int *scratch_i);

SPATIAL_INLINE int spatial_bvh_build_recursive(spatial_bvh *bvh,
                                                int start,
                                                int end,
                                                int depth,
                                                double *scratch_d,
                                                int *scratch_i)
{
    int node_idx = bvh->node_count++;
    
    /* Ensure node array capacity */
    if (node_idx >= bvh->node_capacity) {
        int new_cap = bvh->node_capacity * 2;
        if (new_cap < 64) new_cap = 64;
        spatial_bvh_node *new_nodes = (spatial_bvh_node*)bvh->alloc.malloc(
            sizeof(spatial_bvh_node) * (size_t)new_cap, bvh->alloc.udata);
        if (SPATIAL_UNLIKELY(!new_nodes)) return -1;
        
        if (bvh->nodes) {
            memcpy(new_nodes, bvh->nodes, sizeof(spatial_bvh_node) * (size_t)bvh->node_capacity);
            bvh->alloc.free(bvh->nodes, bvh->alloc.udata);
        }
        bvh->nodes = new_nodes;
        bvh->node_capacity = new_cap;
    }
    
    spatial_bvh_node *node = &bvh->nodes[node_idx];
    spatial_bvh_node_init(node);
    spatial_bvh_compute_bounds(bvh, start, end, node->min, node->max);
    
    int count = end - start;
    
    /* Make leaf if small enough or at max depth */
    if (count <= bvh->leaf_size || depth >= bvh->max_depth) {
        spatial_bvh_node_set_leaf(node, true);
        node->prim_offset = start;
        node->prim_count = count;
        return node_idx;
    }
    
    /* Find best split using SAH */
    spatial_num_t best_cost = SPATIAL_INFINITY;
    int best_axis = -1;
    int best_split = -1;
    
    spatial_num_t total_sa = spatial_bvh_compute_sa(node->min, node->max);
    
    /* Sort by each axis using pre-allocated scratch arrays */
    for (int axis = 0; axis < SPATIAL_BVH_DIMS; axis++) {
        for (int i = 0; i < count; i++) {
            scratch_d[i] = (double)bvh->prims[bvh->indices[start + i]].center[axis];
            scratch_i[i] = bvh->indices[start + i];
        }

        int mid = count / 2;
        spatial_quickselect_double(scratch_d, scratch_i, 0, count - 1, mid);

        /* Copy back to indices for bucket evaluation */
        for (int i = 0; i < count; i++) {
            bvh->indices[start + i] = scratch_i[i];
        }

        /* Build buckets */
        spatial_bvh_bucket buckets[SPATIAL_BVH_NUM_BUCKETS];
        memset(buckets, 0, sizeof(buckets));

        spatial_num_t axis_min = node->min[axis];
        spatial_num_t axis_max = node->max[axis];
        spatial_num_t axis_range = axis_max - axis_min;

        if (axis_range < SPATIAL_EPSILON) continue;

        for (int i = start; i < end; i++) {
            spatial_bvh_prim *p = &bvh->prims[bvh->indices[i]];
            int bucket_idx = (int)((p->center[axis] - axis_min) / axis_range *
                                   (SPATIAL_BVH_NUM_BUCKETS - 0.001));
            if (bucket_idx < 0) bucket_idx = 0;
            if (bucket_idx >= SPATIAL_BVH_NUM_BUCKETS) bucket_idx = SPATIAL_BVH_NUM_BUCKETS - 1;

            spatial_bvh_bucket *b = &buckets[bucket_idx];
            if (b->count == 0) {
                spatial_bbox_copy(p->min, p->max, b->min, b->max, SPATIAL_BVH_DIMS);
            } else {
                for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                    b->min[d] = spatial_min(b->min[d], p->min[d]);
                    b->max[d] = spatial_max(b->max[d], p->max[d]);
                }
            }
            b->count++;
        }

        spatial_num_t cost = spatial_bvh_eval_sah(buckets, SPATIAL_BVH_NUM_BUCKETS, total_sa);

        if (cost < best_cost) {
            best_cost = cost;
            best_axis = axis;
            best_split = start + mid;

            if (best_split <= start || best_split >= end) {
                best_split = start + count / 2;
            }
        }
    }

    if (best_axis < 0 || best_split <= start || best_split >= end) {
        /* Fallback to middle split on longest axis */
        best_split = start + count / 2;
        spatial_num_t best_range = (spatial_num_t)0.0;
        int long_axis = 0;
        for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
            spatial_num_t r = node->max[d] - node->min[d];
            if (r > best_range) {
                best_range = r;
                long_axis = d;
            }
        }
        for (int i = 0; i < count; i++) {
            scratch_d[i] = (double)bvh->prims[bvh->indices[start + i]].center[long_axis];
            scratch_i[i] = bvh->indices[start + i];
        }
        spatial_quickselect_double(scratch_d, scratch_i, 0, count - 1, count / 2);
        for (int i = 0; i < count; i++) {
            bvh->indices[start + i] = scratch_i[i];
        }
    } else {
        /* Re-sort by best axis */
        for (int i = 0; i < count; i++) {
            scratch_d[i] = (double)bvh->prims[bvh->indices[start + i]].center[best_axis];
            scratch_i[i] = bvh->indices[start + i];
        }
        spatial_quickselect_double(scratch_d, scratch_i, 0, count - 1, best_split - start);
        for (int i = 0; i < count; i++) {
            bvh->indices[start + i] = scratch_i[i];
        }
    }

    /* Create children */
    node->left_child  = spatial_bvh_build_recursive(bvh, start,      best_split, depth + 1, scratch_d, scratch_i);
    node->right_child = spatial_bvh_build_recursive(bvh, best_split, end,        depth + 1, scratch_d, scratch_i);
    
    return node_idx;
}

SPATIAL_INLINE bool spatial_bvh_build(spatial_bvh *bvh)
{
    if (SPATIAL_UNLIKELY(!bvh || bvh->prim_count == 0)) return false;
    if (bvh->built) return true;
    
    /* Create indices array */
    bvh->indices = (int*)bvh->alloc.malloc(sizeof(int) * (size_t)bvh->prim_count, bvh->alloc.udata);
    if (SPATIAL_UNLIKELY(!bvh->indices)) return false;
    
    for (int i = 0; i < bvh->prim_count; i++) {
        bvh->indices[i] = i;
    }
    
    /* Allocate initial nodes */
    bvh->node_capacity = bvh->prim_count * 2;
    bvh->nodes = (spatial_bvh_node*)bvh->alloc.malloc(
        sizeof(spatial_bvh_node) * (size_t)bvh->node_capacity, bvh->alloc.udata);
    if (SPATIAL_UNLIKELY(!bvh->nodes)) {
        bvh->alloc.free(bvh->indices, bvh->alloc.udata);
        bvh->indices = NULL;
        return false;
    }
    
    /* Scratch arrays reused across all recursive calls — eliminates O(N log N) mallocs */
    double *scratch_d = (double*)bvh->alloc.malloc(
        sizeof(double) * (size_t)bvh->prim_count, bvh->alloc.udata);
    int *scratch_i = (int*)bvh->alloc.malloc(
        sizeof(int) * (size_t)bvh->prim_count, bvh->alloc.udata);
    if (SPATIAL_UNLIKELY(!scratch_d || !scratch_i)) {
        bvh->alloc.free(scratch_d, bvh->alloc.udata);
        bvh->alloc.free(scratch_i, bvh->alloc.udata);
        bvh->alloc.free(bvh->nodes, bvh->alloc.udata);
        bvh->nodes = NULL;
        bvh->alloc.free(bvh->indices, bvh->alloc.udata);
        bvh->indices = NULL;
        return false;
    }

    bvh->node_count = 0;
    spatial_bvh_build_recursive(bvh, 0, bvh->prim_count, 0, scratch_d, scratch_i);

    bvh->alloc.free(scratch_d, bvh->alloc.udata);
    bvh->alloc.free(scratch_i, bvh->alloc.udata);

    bvh->built = true;
    return true;
}

/* Search */
SPATIAL_INLINE void spatial_bvh_search(const spatial_bvh *bvh,
                                        const spatial_num_t *SPATIAL_RESTRICT qmin,
                                        const spatial_num_t *SPATIAL_RESTRICT qmax,
                                        spatial_bvh_iter_fn iter,
                                        void *udata)
{
    if (SPATIAL_UNLIKELY(!bvh || !bvh->built || !bvh->nodes || !iter)) return;

    int stk[SPATIAL_BVH_MAX_DEPTH * 2 + 4];
    int sp = 0;
    stk[sp++] = 0; /* root */

    while (sp > 0) {
        int node_idx = stk[--sp];
        spatial_bvh_node *node = &bvh->nodes[node_idx];

        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(node->min, node->max, qmin, qmax,
                                                     SPATIAL_BVH_DIMS))) {
            continue;
        }

        if (SPATIAL_LIKELY(spatial_bvh_node_is_leaf(node))) {
            for (int i = 0; i < node->prim_count; i++) {
                int prim_idx = bvh->indices[node->prim_offset + i];
                spatial_bvh_prim *p = &bvh->prims[prim_idx];
                if (spatial_bbox_overlaps(p->min, p->max, qmin, qmax, SPATIAL_BVH_DIMS)) {
                    if (!iter(p->min, p->max, p->data, udata)) return;
                }
            }
        } else {
            if (node->left_child  >= 0) stk[sp++] = node->left_child;
            if (node->right_child >= 0) stk[sp++] = node->right_child;
        }
    }
}

/* Scan all items */
SPATIAL_INLINE void spatial_bvh_scan(const spatial_bvh *bvh,
                                      spatial_bvh_iter_fn iter,
                                      void *udata)
{
    if (SPATIAL_UNLIKELY(!bvh || !bvh->built || !bvh->nodes || !iter)) return;

    int stk[SPATIAL_BVH_MAX_DEPTH * 2 + 4];
    int sp = 0;
    stk[sp++] = 0;

    while (sp > 0) {
        int node_idx = stk[--sp];
        spatial_bvh_node *node = &bvh->nodes[node_idx];

        if (SPATIAL_LIKELY(spatial_bvh_node_is_leaf(node))) {
            for (int i = 0; i < node->prim_count; i++) {
                int prim_idx = bvh->indices[node->prim_offset + i];
                spatial_bvh_prim *p = &bvh->prims[prim_idx];
                if (!iter(p->min, p->max, p->data, udata)) return;
            }
        } else {
            if (node->left_child  >= 0) stk[sp++] = node->left_child;
            if (node->right_child >= 0) stk[sp++] = node->right_child;
        }
    }
}

SPATIAL_INLINE int spatial_bvh_count(const spatial_bvh *bvh) {
    return bvh ? bvh->prim_count : 0;
}

/* Ray intersection */
SPATIAL_INLINE bool spatial_bvh_intersect_ray(const spatial_bvh *bvh,
                                               const spatial_bvh_ray *ray,
                                               spatial_bvh_hit *out_hit)
{
    if (SPATIAL_UNLIKELY(!bvh || !bvh->built || !bvh->nodes || !out_hit)) return false;

    spatial_num_t closest_t = ray->tmax;
    bool hit = false;

    int stk[SPATIAL_BVH_MAX_DEPTH * 2 + 4];
    int sp = 0;
    stk[sp++] = 0;
    
    while (sp > 0) {
        int node_idx = stk[--sp];
        spatial_bvh_node *node = &bvh->nodes[node_idx];

        if (!spatial_bvh_ray_intersect_box(ray, node->min, node->max, ray->tmin, closest_t)) {
            continue;
        }

        if (spatial_bvh_node_is_leaf(node)) {
            for (int i = 0; i < node->prim_count; i++) {
                int prim_idx = bvh->indices[node->prim_offset + i];
                spatial_bvh_prim *p = &bvh->prims[prim_idx];

                if (spatial_bvh_ray_intersect_box(ray, p->min, p->max, ray->tmin, closest_t)) {
                    spatial_num_t t_entry = (spatial_num_t)0.0;
                    for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                        spatial_num_t t1 = (p->min[d] - ray->origin[d]) * ray->inv_dir[d];
                        spatial_num_t t2 = (p->max[d] - ray->origin[d]) * ray->inv_dir[d];
                        spatial_num_t tn = spatial_min(t1, t2);
                        t_entry = spatial_max(t_entry, tn);
                    }

                    if (t_entry >= ray->tmin && t_entry < closest_t) {
                        closest_t = t_entry;
                        out_hit->data = p->data;
                        out_hit->t = t_entry;
                        out_hit->prim_id = prim_idx;
                        hit = true;
                    }
                }
            }
        } else {
            /* Push children in far-to-near order for better culling */
            spatial_num_t left_t = closest_t;
            spatial_num_t right_t = closest_t;
            bool left_hit = false, right_hit = false;

            if (node->left_child >= 0) {
                spatial_bvh_node *left = &bvh->nodes[node->left_child];
                if (spatial_bvh_ray_intersect_box(ray, left->min, left->max, ray->tmin, closest_t)) {
                    spatial_num_t center[SPATIAL_BVH_DIMS];
                    for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                        center[d] = (left->min[d] + left->max[d]) * (spatial_num_t)0.5;
                    }
                    left_t = (spatial_num_t)0.0;
                    for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                        spatial_num_t diff = center[d] - ray->origin[d];
                        left_t += diff * diff;
                    }
                    left_hit = true;
                }
            }

            if (node->right_child >= 0) {
                spatial_bvh_node *right = &bvh->nodes[node->right_child];
                if (spatial_bvh_ray_intersect_box(ray, right->min, right->max, ray->tmin, closest_t)) {
                    spatial_num_t center[SPATIAL_BVH_DIMS];
                    for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                        center[d] = (right->min[d] + right->max[d]) * (spatial_num_t)0.5;
                    }
                    right_t = (spatial_num_t)0.0;
                    for (int d = 0; d < SPATIAL_BVH_DIMS; d++) {
                        spatial_num_t diff = center[d] - ray->origin[d];
                        right_t += diff * diff;
                    }
                    right_hit = true;
                }
            }

            if (left_hit && right_hit) {
                if (left_t > right_t) {
                    stk[sp++] = node->left_child;
                    stk[sp++] = node->right_child;
                } else {
                    stk[sp++] = node->right_child;
                    stk[sp++] = node->left_child;
                }
            } else if (left_hit) {
                stk[sp++] = node->left_child;
            } else if (right_hit) {
                stk[sp++] = node->right_child;
            }
        }
    }

    return hit;
}

/* Delete (not efficient for BVH - requires rebuild) */
SPATIAL_INLINE bool spatial_bvh_delete(spatial_bvh *bvh,
                                        const spatial_num_t *min,
                                        const spatial_num_t *max,
                                        spatial_data_t data)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!bvh || bvh->built)) return false;
    
    for (int i = 0; i < bvh->prim_count; i++) {
        if (bvh->prims[i].data == data) {
            if (bvh->callbacks.free) {
                bvh->callbacks.free(bvh->prims[i].data, bvh->udata);
            }
            bvh->prims[i] = bvh->prims[bvh->prim_count - 1];
            bvh->prim_count--;
            return true;
        }
    }
    
    return false;
}

typedef bool (*spatial_bvh_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_bvh_delete_with_comparator(spatial_bvh *bvh,
                                                        const spatial_num_t *min,
                                                        const spatial_num_t *max,
                                                        spatial_bvh_cmp_fn cmp,
                                                        void *cmp_udata)
{
    (void)min; (void)max;
    if (SPATIAL_UNLIKELY(!bvh || bvh->built)) return false;
    
    for (int i = 0; i < bvh->prim_count; i++) {
        if (cmp(bvh->prims[i].data, cmp_udata)) {
            if (bvh->callbacks.free) {
                bvh->callbacks.free(bvh->prims[i].data, bvh->udata);
            }
            bvh->prims[i] = bvh->prims[bvh->prim_count - 1];
            bvh->prim_count--;
            return true;
        }
    }
    
    return false;
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_BVH_H */
