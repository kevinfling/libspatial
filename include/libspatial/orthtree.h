/**
 * @file   orthtree.h
 * @brief  N-Dimensional Orthtree for spatial indexing
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * A orthtree recursively subdivides N-dimensional space into 2^N orthants.
 * Configurable dimensionality at compile time via SPATIAL_ORTHTREE_DIMS.
 */

#ifndef SPATIAL_ORTHTREE_H
#define SPATIAL_ORTHTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

#ifndef SPATIAL_ORTHTREE_DIMS
    #define SPATIAL_ORTHTREE_DIMS 4
#endif

#define SPATIAL_ORTHTREE_CHILDREN (1 << SPATIAL_ORTHTREE_DIMS)

#ifndef SPATIAL_ORTHTREE_MAX_DEPTH
    #define SPATIAL_ORTHTREE_MAX_DEPTH 12
#endif

#ifndef SPATIAL_ORTHTREE_MAX_ITEMS
    #define SPATIAL_ORTHTREE_MAX_ITEMS 8
#endif

SPATIAL_STATIC_ASSERT(orthtree_dims_valid,
    SPATIAL_ORTHTREE_DIMS > 0 && SPATIAL_ORTHTREE_DIMS <= 8);
SPATIAL_STATIC_ASSERT(orthtree_children_valid,
    SPATIAL_ORTHTREE_CHILDREN > 0 && SPATIAL_ORTHTREE_CHILDREN <= 256);

/* Forward declarations */
struct spatial_orthtree_node;

/* Leaf storage: Structure-of-Arrays with exponential growth */
typedef struct spatial_orthtree_leaf {
    spatial_num_t *mins;      /* capacity * dims */
    spatial_num_t *maxs;      /* capacity * dims */
    spatial_data_t *datas;    /* capacity */
    int count;
    int capacity;
} spatial_orthtree_leaf;

/* Node structure (cache-line aligned for SIMD-ready traversal) */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_orthtree_node {
    spatial_num_t min[SPATIAL_ORTHTREE_DIMS];
    spatial_num_t max[SPATIAL_ORTHTREE_DIMS];
    spatial_num_t center[SPATIAL_ORTHTREE_DIMS];
    
    union {
        struct spatial_orthtree_node *children[SPATIAL_ORTHTREE_CHILDREN];
        spatial_orthtree_leaf leaf;
    };
    
    int depth;
    bool is_leaf;
} spatial_orthtree_node;

/* Main orthtree structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_orthtree {
    spatial_orthtree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    int count;
    int max_depth;
    int max_items;
    bool relaxed_atomics;
} spatial_orthtree;

/* Iterator callback */
typedef bool (*spatial_orthtree_iter_fn)(const spatial_num_t *min,
                                             const spatial_num_t *max,
                                             spatial_data_t data,
                                             void *udata);

/* Forward declarations for public API */
SPATIAL_INLINE spatial_orthtree* spatial_orthtree_new_with_allocator(const spatial_allocator *alloc);

/* ============================================================================
 * Internal Node Management
 * ============================================================================ */

SPATIAL_INLINE spatial_orthtree_node* spatial_orthtree_node_new(
    const spatial_num_t *SPATIAL_RESTRICT min,
    const spatial_num_t *SPATIAL_RESTRICT max,
    int depth,
    const spatial_allocator *alloc)
{
    spatial_orthtree_node *node = (spatial_orthtree_node*)alloc->malloc(
        sizeof(spatial_orthtree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    spatial_bbox_copy(min, max, node->min, node->max, SPATIAL_ORTHTREE_DIMS);
    
    for (int i = 0; i < SPATIAL_ORTHTREE_DIMS; i++) {
        node->center[i] = (min[i] + max[i]) * (spatial_num_t)0.5;
    }
    
    for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
        node->children[i] = NULL;
    }
    
    node->depth = depth;
    node->is_leaf = true;
    node->leaf.mins = NULL;
    node->leaf.maxs = NULL;
    node->leaf.datas = NULL;
    node->leaf.count = 0;
    node->leaf.capacity = 0;
    
    return node;
}

SPATIAL_INLINE void spatial_orthtree_leaf_free(spatial_orthtree_leaf *leaf,
                                                    const spatial_allocator *alloc,
                                                    const spatial_item_callbacks *cb,
                                                    void *udata)
{
    if (cb && cb->free) {
        for (int i = 0; i < leaf->count; i++) {
            cb->free(leaf->datas[i], udata);
        }
    }
    if (leaf->mins)    alloc->free(leaf->mins, alloc->udata);
    if (leaf->maxs)    alloc->free(leaf->maxs, alloc->udata);
    if (leaf->datas)   alloc->free(leaf->datas, alloc->udata);
    leaf->mins = NULL;
    leaf->maxs = NULL;
    leaf->datas = NULL;
    leaf->count = 0;
    leaf->capacity = 0;
}

SPATIAL_INLINE void spatial_orthtree_node_free(spatial_orthtree_node *node,
                                               const spatial_allocator *alloc,
                                               const spatial_item_callbacks *cb,
                                               void *udata)
{
    if (SPATIAL_UNLIKELY(!node)) return;
    
    if (node->is_leaf) {
        spatial_orthtree_leaf_free(&node->leaf, alloc, cb, udata);
    } else {
        for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
            if (node->children[i]) {
                spatial_orthtree_node_free(node->children[i], alloc, cb, udata);
            }
        }
    }
    alloc->free(node, alloc->udata);
}

SPATIAL_INLINE int spatial_orthtree_octant(const spatial_orthtree_node *SPATIAL_RESTRICT node,
                                               const spatial_num_t *SPATIAL_RESTRICT min,
                                               const spatial_num_t *SPATIAL_RESTRICT max)
{
    int octant = 0;
    for (int i = 0; i < SPATIAL_ORTHTREE_DIMS; i++) {
        spatial_num_t pc = (min[i] + max[i]) * (spatial_num_t)0.5;
        if (pc >= node->center[i]) {
            octant |= (1 << i);
        }
    }
    return octant;
}

SPATIAL_INLINE void spatial_orthtree_child_bounds(const spatial_orthtree_node *SPATIAL_RESTRICT node,
                                                      int octant,
                                                      spatial_num_t *SPATIAL_RESTRICT child_min,
                                                      spatial_num_t *SPATIAL_RESTRICT child_max)
{
    for (int i = 0; i < SPATIAL_ORTHTREE_DIMS; i++) {
        if (octant & (1 << i)) {
            child_min[i] = node->center[i];
            child_max[i] = node->max[i];
        } else {
            child_min[i] = node->min[i];
            child_max[i] = node->center[i];
        }
    }
}

SPATIAL_INLINE bool spatial_orthtree_leaf_grow(spatial_orthtree_leaf *leaf,
                                                     int dims,
                                                     const spatial_allocator *alloc)
{
    int new_cap = leaf->capacity ? leaf->capacity * 2 : 4;
    spatial_num_t *nm = (spatial_num_t*)alloc->malloc(
        sizeof(spatial_num_t) * new_cap * dims, alloc->udata);
    spatial_num_t *nx = (spatial_num_t*)alloc->malloc(
        sizeof(spatial_num_t) * new_cap * dims, alloc->udata);
    spatial_data_t *nd = (spatial_data_t*)alloc->malloc(
        sizeof(spatial_data_t) * new_cap, alloc->udata);
    if (SPATIAL_UNLIKELY(!nm || !nx || !nd)) {
        if (nm) alloc->free(nm, alloc->udata);
        if (nx) alloc->free(nx, alloc->udata);
        if (nd) alloc->free(nd, alloc->udata);
        return false;
    }
    if (leaf->count > 0) {
        memcpy(nm, leaf->mins, sizeof(spatial_num_t) * leaf->count * dims);
        memcpy(nx, leaf->maxs, sizeof(spatial_num_t) * leaf->count * dims);
        memcpy(nd, leaf->datas, sizeof(spatial_data_t) * leaf->count);
    }
    if (leaf->mins)  alloc->free(leaf->mins, alloc->udata);
    if (leaf->maxs)  alloc->free(leaf->maxs, alloc->udata);
    if (leaf->datas) alloc->free(leaf->datas, alloc->udata);
    leaf->mins = nm;
    leaf->maxs = nx;
    leaf->datas = nd;
    leaf->capacity = new_cap;
    return true;
}

SPATIAL_INLINE bool spatial_orthtree_node_split(spatial_orthtree_node *SPATIAL_RESTRICT node,
                                                int max_depth,
                                                const spatial_allocator *alloc)
{
    if (SPATIAL_UNLIKELY(node->depth >= max_depth)) return false;
    
    /* Save old leaf data before union is overwritten by children array */
    spatial_orthtree_leaf old_leaf = node->leaf;
    
    for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
        node->children[i] = NULL;
    }
    node->is_leaf = false;
    
    spatial_num_t cmin[SPATIAL_ORTHTREE_DIMS];
    spatial_num_t cmax[SPATIAL_ORTHTREE_DIMS];
    
    for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
        spatial_orthtree_child_bounds(node, i, cmin, cmax);
        node->children[i] = spatial_orthtree_node_new(cmin, cmax, node->depth + 1, alloc);
        if (SPATIAL_UNLIKELY(!node->children[i])) {
            for (int j = 0; j < i; j++) {
                spatial_orthtree_node_free(node->children[j], alloc, NULL, NULL);
                node->children[j] = NULL;
            }
            node->is_leaf = true;
            node->leaf = old_leaf;
            return false;
        }
    }
    
    /* Redistribute items from old leaf arrays into children */
    for (int i = 0; i < old_leaf.count; i++) {
        int oct = spatial_orthtree_octant(node,
            &old_leaf.mins[i * SPATIAL_ORTHTREE_DIMS],
            &old_leaf.maxs[i * SPATIAL_ORTHTREE_DIMS]);
        spatial_orthtree_node *child = node->children[oct];
        spatial_orthtree_leaf *cl = &child->leaf;
        
        if (SPATIAL_UNLIKELY(cl->count >= cl->capacity)) {
            if (SPATIAL_UNLIKELY(!spatial_orthtree_leaf_grow(cl, SPATIAL_ORTHTREE_DIMS, alloc))) {
                /* OOM during redistribution: restore old leaf, free children */
                for (int j = 0; j < SPATIAL_ORTHTREE_CHILDREN; j++) {
                    if (node->children[j]) {
                        spatial_orthtree_node_free(node->children[j], alloc, NULL, NULL);
                        node->children[j] = NULL;
                    }
                }
                node->is_leaf = true;
                node->leaf = old_leaf;
                return false;
            }
        }
        memcpy(&cl->mins[cl->count * SPATIAL_ORTHTREE_DIMS],
               &old_leaf.mins[i * SPATIAL_ORTHTREE_DIMS],
               sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
        memcpy(&cl->maxs[cl->count * SPATIAL_ORTHTREE_DIMS],
               &old_leaf.maxs[i * SPATIAL_ORTHTREE_DIMS],
               sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
        cl->datas[cl->count] = old_leaf.datas[i];
        cl->count++;
    }
    
    /* Free old leaf arrays */
    if (old_leaf.mins)  alloc->free(old_leaf.mins, alloc->udata);
    if (old_leaf.maxs)  alloc->free(old_leaf.maxs, alloc->udata);
    if (old_leaf.datas) alloc->free(old_leaf.datas, alloc->udata);
    
    return true;
}

/* ============================================================================
 * Public API
 * ============================================================================ */

SPATIAL_INLINE spatial_orthtree* spatial_orthtree_new(void) {
    return spatial_orthtree_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_orthtree* spatial_orthtree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();
    
    spatial_orthtree *ht = (spatial_orthtree*)alloc->malloc(sizeof(spatial_orthtree), alloc->udata);
    if (SPATIAL_UNLIKELY(!ht)) return NULL;
    
    memset(ht, 0, sizeof(spatial_orthtree));
    spatial_refcount_init(&ht->refcount);
    ht->alloc = *alloc;
    ht->max_depth = SPATIAL_ORTHTREE_MAX_DEPTH;
    ht->max_items = SPATIAL_ORTHTREE_MAX_ITEMS;
    ht->relaxed_atomics = true;
    ht->root = NULL;
    
    return ht;
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_orthtree* spatial_orthtree_new_with_memento_arena(memento_arena_t *arena) {
    spatial_allocator alloc = spatial_memento_arena_allocator(arena);
    return spatial_orthtree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_orthtree* spatial_orthtree_new_with_memento_pool(memento_pool_t *pool) {
    spatial_allocator alloc = spatial_memento_pool_allocator(pool);
    return spatial_orthtree_new_with_allocator(&alloc);
}
#endif

SPATIAL_INLINE spatial_orthtree* spatial_orthtree_new_with_arena(spatial_arena_t *arena) {
    spatial_allocator alloc = spatial_arena_allocator(arena);
    return spatial_orthtree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_orthtree* spatial_orthtree_new_with_pool(spatial_pool_t *pool) {
    spatial_allocator alloc = spatial_pool_allocator(pool);
    return spatial_orthtree_new_with_allocator(&alloc);
}

SPATIAL_INLINE void spatial_orthtree_free(spatial_orthtree *ht) {
    if (SPATIAL_UNLIKELY(!ht)) return;
    
    if (spatial_refcount_decrement(&ht->refcount) > 0) return;
    
    spatial_orthtree_node_free(ht->root, &ht->alloc, &ht->callbacks, ht->udata);
    ht->alloc.free(ht, ht->alloc.udata);
}

SPATIAL_INLINE spatial_orthtree* spatial_orthtree_clone(const spatial_orthtree *ht) {
    if (SPATIAL_UNLIKELY(!ht)) return NULL;
    
    spatial_refcount_increment(&ht->refcount);
    union { const spatial_orthtree *in; spatial_orthtree *out; } pun; pun.in = ht; return pun.out;
}

SPATIAL_INLINE void spatial_orthtree_set_item_callbacks(spatial_orthtree *ht,
                                                            const spatial_item_callbacks *cb)
{
    if (ht) ht->callbacks = *cb;
}

SPATIAL_INLINE void spatial_orthtree_set_udata(spatial_orthtree *ht, void *udata) {
    if (ht) ht->udata = udata;
}

SPATIAL_INLINE void spatial_orthtree_opt_relaxed_atomics(spatial_orthtree *ht, bool enable) {
    if (ht) ht->relaxed_atomics = enable;
}

SPATIAL_INLINE bool spatial_orthtree_init_root(spatial_orthtree *ht,
                                                   const spatial_num_t *min,
                                                   const spatial_num_t *max)
{
    if (ht->root) return true;
    
    ht->root = spatial_orthtree_node_new(min, max, 0, &ht->alloc);
    return ht->root != NULL;
}


SPATIAL_INLINE bool spatial_orthtree_expand_root(spatial_orthtree *ht,
                                                     const spatial_num_t *min,
                                                     const spatial_num_t *max)
{
    if (!ht->root) return spatial_orthtree_init_root(ht, min, max);
    
    spatial_num_t new_min[SPATIAL_ORTHTREE_DIMS];
    spatial_num_t new_max[SPATIAL_ORTHTREE_DIMS];
    
    for (int i = 0; i < SPATIAL_ORTHTREE_DIMS; i++) {
        new_min[i] = spatial_min(ht->root->min[i], min[i]);
        new_max[i] = spatial_max(ht->root->max[i], max[i]);
        spatial_num_t size = new_max[i] - new_min[i];
        if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
        new_min[i] -= size * (spatial_num_t)0.05;
        new_max[i] += size * (spatial_num_t)0.05;
    }
    
    if (ht->root->is_leaf) {
        for (int i = 0; i < SPATIAL_ORTHTREE_DIMS; i++) {
            ht->root->min[i] = new_min[i];
            ht->root->max[i] = new_max[i];
            ht->root->center[i] = (new_min[i] + new_max[i]) * (spatial_num_t)0.5;
        }
        return true;
    }
    
    /* Root has children: create new parent root */
    spatial_orthtree_node *new_root = spatial_orthtree_node_new(new_min, new_max, 0, &ht->alloc);
    if (!new_root) return false;
    
    new_root->is_leaf = false;
    for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
        new_root->children[i] = NULL;
    }
    
    int old_octant = spatial_orthtree_octant(new_root, ht->root->min, ht->root->max);
    new_root->children[old_octant] = ht->root;
    ht->root = new_root;
    return true;
}

SPATIAL_INLINE bool spatial_orthtree_insert(spatial_orthtree *ht,
                                                const spatial_num_t *SPATIAL_RESTRICT min,
                                                const spatial_num_t *SPATIAL_RESTRICT max,
                                                spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!ht)) return false;
    
    if (SPATIAL_UNLIKELY(!ht->root)) {
        spatial_num_t root_min[SPATIAL_ORTHTREE_DIMS];
        spatial_num_t root_max[SPATIAL_ORTHTREE_DIMS];
        
        for (int i = 0; i < SPATIAL_ORTHTREE_DIMS; i++) {
            root_min[i] = min[i];
            root_max[i] = max[i];
        }
        
        for (int i = 0; i < SPATIAL_ORTHTREE_DIMS; i++) {
            spatial_num_t size = root_max[i] - root_min[i];
            if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
            root_min[i] -= size * (spatial_num_t)0.1;
            root_max[i] += size * (spatial_num_t)0.1;
        }
        
        bool ok = spatial_orthtree_init_root(ht, root_min, root_max);
        if (!ok) return false;
    }
    
    if (!spatial_bbox_contains(ht->root->min, ht->root->max, min, max, SPATIAL_ORTHTREE_DIMS)) {
        if (!spatial_orthtree_expand_root(ht, min, max)) return false;
    }
    
    spatial_orthtree_node *node = ht->root;
    
    while (!node->is_leaf) {
        int octant = spatial_orthtree_octant(node, min, max);
        if (!node->children[octant]) {
            spatial_num_t child_min[SPATIAL_ORTHTREE_DIMS];
            spatial_num_t child_max[SPATIAL_ORTHTREE_DIMS];
            spatial_orthtree_child_bounds(node, octant, child_min, child_max);
            node->children[octant] = spatial_orthtree_node_new(child_min, child_max, node->depth + 1, &ht->alloc);
            if (SPATIAL_UNLIKELY(!node->children[octant])) return false;
        }
        node = node->children[octant];
    }
    
    spatial_orthtree_leaf *leaf = &node->leaf;
    if (SPATIAL_UNLIKELY(leaf->count >= leaf->capacity)) {
        if (SPATIAL_UNLIKELY(!spatial_orthtree_leaf_grow(leaf, SPATIAL_ORTHTREE_DIMS, &ht->alloc))) {
            return false;
        }
    }
    memcpy(&leaf->mins[leaf->count * SPATIAL_ORTHTREE_DIMS],
           min, sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
    memcpy(&leaf->maxs[leaf->count * SPATIAL_ORTHTREE_DIMS],
           max, sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
    leaf->datas[leaf->count] = data;
    leaf->count++;
    ht->count++;
    
    if (leaf->count > ht->max_items && node->depth < ht->max_depth) {
        spatial_orthtree_node_split(node, ht->max_depth, &ht->alloc);
    }
    
    return true;
}

SPATIAL_INLINE void spatial_orthtree_search(const spatial_orthtree *ht,
                                                const spatial_num_t *SPATIAL_RESTRICT qmin,
                                                const spatial_num_t *SPATIAL_RESTRICT qmax,
                                                spatial_orthtree_iter_fn iter,
                                                void *udata)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root || !iter)) return;

    spatial_orthtree_node *stk[(SPATIAL_ORTHTREE_MAX_DEPTH + 1) * SPATIAL_ORTHTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = ht->root;

    while (sp > 0) {
        spatial_orthtree_node *node = stk[--sp];

        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(node->min, node->max, qmin, qmax, SPATIAL_ORTHTREE_DIMS))) {
            continue;
        }

        if (node->is_leaf) {
            const spatial_orthtree_leaf *leaf = &node->leaf;
            for (int i = 0; i < leaf->count; i++) {
                if (spatial_bbox_overlaps(&leaf->mins[i * SPATIAL_ORTHTREE_DIMS],
                                          &leaf->maxs[i * SPATIAL_ORTHTREE_DIMS],
                                          qmin, qmax, SPATIAL_ORTHTREE_DIMS)) {
                    if (!iter(&leaf->mins[i * SPATIAL_ORTHTREE_DIMS],
                              &leaf->maxs[i * SPATIAL_ORTHTREE_DIMS],
                              leaf->datas[i], udata)) {
                        return;
                    }
                }
            }
        } else {
            for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE void spatial_orthtree_scan(const spatial_orthtree *ht,
                                              spatial_orthtree_iter_fn iter,
                                              void *udata)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root || !iter)) return;

    spatial_orthtree_node *stk[(SPATIAL_ORTHTREE_MAX_DEPTH + 1) * SPATIAL_ORTHTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = ht->root;

    while (sp > 0) {
        spatial_orthtree_node *node = stk[--sp];

        if (node->is_leaf) {
            const spatial_orthtree_leaf *leaf = &node->leaf;
            for (int i = 0; i < leaf->count; i++) {
                if (!iter(&leaf->mins[i * SPATIAL_ORTHTREE_DIMS],
                          &leaf->maxs[i * SPATIAL_ORTHTREE_DIMS],
                          leaf->datas[i], udata)) {
                    return;
                }
            }
        } else {
            for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE int spatial_orthtree_count(const spatial_orthtree *ht) {
    return ht ? ht->count : 0;
}

SPATIAL_INLINE bool spatial_orthtree_delete(spatial_orthtree *ht,
                                                const spatial_num_t *min,
                                                const spatial_num_t *max,
                                                spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root)) return false;

    spatial_orthtree_node *stk[(SPATIAL_ORTHTREE_MAX_DEPTH + 1) * SPATIAL_ORTHTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = ht->root;

    while (sp > 0 && !found) {
        spatial_orthtree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_ORTHTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_orthtree_leaf *leaf = &node->leaf;
            for (int i = 0; i < leaf->count; i++) {
                if (leaf->datas[i] == data &&
                    spatial_bbox_overlaps(&leaf->mins[i * SPATIAL_ORTHTREE_DIMS],
                                          &leaf->maxs[i * SPATIAL_ORTHTREE_DIMS],
                                          min, max, SPATIAL_ORTHTREE_DIMS)) {
                    if (ht->callbacks.free) {
                        ht->callbacks.free(leaf->datas[i], ht->udata);
                    }
                    int last = leaf->count - 1;
                    if (i != last) {
                        memcpy(&leaf->mins[i * SPATIAL_ORTHTREE_DIMS],
                               &leaf->mins[last * SPATIAL_ORTHTREE_DIMS],
                               sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
                        memcpy(&leaf->maxs[i * SPATIAL_ORTHTREE_DIMS],
                               &leaf->maxs[last * SPATIAL_ORTHTREE_DIMS],
                               sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
                        leaf->datas[i] = leaf->datas[last];
                    }
                    leaf->count--;
                    ht->count--;
                    found = true;
                    break;
                }
            }
        } else {
            for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }

    return found;
}

typedef bool (*spatial_orthtree_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_orthtree_delete_with_comparator(spatial_orthtree *ht,
                                                            const spatial_num_t *min,
                                                            const spatial_num_t *max,
                                                            spatial_orthtree_cmp_fn cmp,
                                                            void *cmp_udata)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root)) return false;

    spatial_orthtree_node *stk[(SPATIAL_ORTHTREE_MAX_DEPTH + 1) * SPATIAL_ORTHTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = ht->root;

    while (sp > 0 && !found) {
        spatial_orthtree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_ORTHTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_orthtree_leaf *leaf = &node->leaf;
            for (int i = 0; i < leaf->count; i++) {
                if (cmp(leaf->datas[i], cmp_udata) &&
                    spatial_bbox_overlaps(&leaf->mins[i * SPATIAL_ORTHTREE_DIMS],
                                          &leaf->maxs[i * SPATIAL_ORTHTREE_DIMS],
                                          min, max, SPATIAL_ORTHTREE_DIMS)) {
                    if (ht->callbacks.free) {
                        ht->callbacks.free(leaf->datas[i], ht->udata);
                    }
                    int last = leaf->count - 1;
                    if (i != last) {
                        memcpy(&leaf->mins[i * SPATIAL_ORTHTREE_DIMS],
                               &leaf->mins[last * SPATIAL_ORTHTREE_DIMS],
                               sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
                        memcpy(&leaf->maxs[i * SPATIAL_ORTHTREE_DIMS],
                               &leaf->maxs[last * SPATIAL_ORTHTREE_DIMS],
                               sizeof(spatial_num_t) * SPATIAL_ORTHTREE_DIMS);
                        leaf->datas[i] = leaf->datas[last];
                    }
                    leaf->count--;
                    ht->count--;
                    found = true;
                    break;
                }
            }
        } else {
            for (int i = 0; i < SPATIAL_ORTHTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }

    return found;
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_ORTHTREE_H */
