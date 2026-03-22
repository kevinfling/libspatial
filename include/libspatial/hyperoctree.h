/**
 * @file   hyperoctree.h
 * @brief  N-Dimensional Hyperoctree for spatial indexing
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * A hyperoctree recursively subdivides N-dimensional space into 2^N hyperoctants.
 * Configurable dimensionality at compile time via SPATIAL_HYPEROCTREE_DIMS.
 */

#ifndef SPATIAL_HYPEROCTREE_H
#define SPATIAL_HYPEROCTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

#ifndef SPATIAL_HYPEROCTREE_DIMS
    #define SPATIAL_HYPEROCTREE_DIMS 4
#endif

#define SPATIAL_HYPEROCTREE_CHILDREN (1 << SPATIAL_HYPEROCTREE_DIMS)

#ifndef SPATIAL_HYPEROCTREE_MAX_DEPTH
    #define SPATIAL_HYPEROCTREE_MAX_DEPTH 12
#endif

#ifndef SPATIAL_HYPEROCTREE_MAX_ITEMS
    #define SPATIAL_HYPEROCTREE_MAX_ITEMS 8
#endif

SPATIAL_STATIC_ASSERT(hyperoctree_dims_valid,
    SPATIAL_HYPEROCTREE_DIMS > 0 && SPATIAL_HYPEROCTREE_DIMS <= 8);
SPATIAL_STATIC_ASSERT(hyperoctree_children_valid,
    SPATIAL_HYPEROCTREE_CHILDREN > 0 && SPATIAL_HYPEROCTREE_CHILDREN <= 256);

/* Forward declarations */
struct spatial_hyperoctree_node;

/* Item structure - fixed size, no extra malloc */
typedef struct spatial_hyperoctree_item {
    spatial_num_t min[SPATIAL_HYPEROCTREE_DIMS];
    spatial_num_t max[SPATIAL_HYPEROCTREE_DIMS];
    spatial_data_t data;
    struct spatial_hyperoctree_item *next;
} spatial_hyperoctree_item;

/* Node structure (cache-line aligned for SIMD-ready traversal) */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_hyperoctree_node {
    spatial_num_t min[SPATIAL_HYPEROCTREE_DIMS];
    spatial_num_t max[SPATIAL_HYPEROCTREE_DIMS];
    spatial_num_t center[SPATIAL_HYPEROCTREE_DIMS];
    
    union {
        struct spatial_hyperoctree_node *children[SPATIAL_HYPEROCTREE_CHILDREN];
        spatial_hyperoctree_item *items;
    };
    
    int item_count;
    int depth;
    bool is_leaf;
} spatial_hyperoctree_node;

/* Main hyperoctree structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_hyperoctree {
    spatial_hyperoctree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    int count;
    int max_depth;
    int max_items;
    bool relaxed_atomics;
} spatial_hyperoctree;

/* Iterator callback */
typedef bool (*spatial_hyperoctree_iter_fn)(const spatial_num_t *min,
                                             const spatial_num_t *max,
                                             spatial_data_t data,
                                             void *udata);

/* Forward declarations for public API */
SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_new_with_allocator(const spatial_allocator *alloc);

/* ============================================================================
 * Internal Node Management
 * ============================================================================ */

SPATIAL_INLINE spatial_hyperoctree_node* spatial_hyperoctree_node_new(
    const spatial_num_t *SPATIAL_RESTRICT min,
    const spatial_num_t *SPATIAL_RESTRICT max,
    int depth,
    const spatial_allocator *alloc)
{
    spatial_hyperoctree_node *node = (spatial_hyperoctree_node*)alloc->malloc(
        sizeof(spatial_hyperoctree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    spatial_bbox_copy(min, max, node->min, node->max, SPATIAL_HYPEROCTREE_DIMS);
    
    for (int i = 0; i < SPATIAL_HYPEROCTREE_DIMS; i++) {
        node->center[i] = (min[i] + max[i]) * (spatial_num_t)0.5;
    }
    
    for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
        node->children[i] = NULL;
    }
    
    node->item_count = 0;
    node->depth = depth;
    node->is_leaf = true;
    node->items = NULL;
    
    return node;
}

SPATIAL_INLINE void spatial_hyperoctree_node_free(spatial_hyperoctree_node *node,
                                                   const spatial_allocator *alloc,
                                                   const spatial_item_callbacks *cb,
                                                   void *udata);

SPATIAL_INLINE void spatial_hyperoctree_item_free_list(spatial_hyperoctree_item *item,
                                                        const spatial_allocator *alloc,
                                                        const spatial_item_callbacks *cb,
                                                        void *udata)
{
    while (item) {
        spatial_hyperoctree_item *next = item->next;
        if (cb && cb->free) {
            cb->free(item->data, udata);
        }
        alloc->free(item, alloc->udata);
        item = next;
    }
}

SPATIAL_INLINE void spatial_hyperoctree_node_free(spatial_hyperoctree_node *node,
                                                   const spatial_allocator *alloc,
                                                   const spatial_item_callbacks *cb,
                                                   void *udata)
{
    if (SPATIAL_UNLIKELY(!node)) return;
    
    if (node->is_leaf) {
        spatial_hyperoctree_item_free_list(node->items, alloc, cb, udata);
    } else {
        for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
            if (node->children[i]) {
                spatial_hyperoctree_node_free(node->children[i], alloc, cb, udata);
            }
        }
    }
    alloc->free(node, alloc->udata);
}

SPATIAL_INLINE int spatial_hyperoctree_octant(const spatial_hyperoctree_node *SPATIAL_RESTRICT node,
                                               const spatial_num_t *SPATIAL_RESTRICT min,
                                               const spatial_num_t *SPATIAL_RESTRICT max)
{
    int octant = 0;
    for (int i = 0; i < SPATIAL_HYPEROCTREE_DIMS; i++) {
        spatial_num_t pc = (min[i] + max[i]) * (spatial_num_t)0.5;
        if (pc >= node->center[i]) {
            octant |= (1 << i);
        }
    }
    return octant;
}

SPATIAL_INLINE void spatial_hyperoctree_child_bounds(const spatial_hyperoctree_node *SPATIAL_RESTRICT node,
                                                      int octant,
                                                      spatial_num_t *SPATIAL_RESTRICT child_min,
                                                      spatial_num_t *SPATIAL_RESTRICT child_max)
{
    for (int i = 0; i < SPATIAL_HYPEROCTREE_DIMS; i++) {
        if (octant & (1 << i)) {
            child_min[i] = node->center[i];
            child_max[i] = node->max[i];
        } else {
            child_min[i] = node->min[i];
            child_max[i] = node->center[i];
        }
    }
}

SPATIAL_INLINE bool spatial_hyperoctree_node_split(spatial_hyperoctree_node *SPATIAL_RESTRICT node,
                                                    int max_depth,
                                                    const spatial_allocator *alloc)
{
    if (SPATIAL_UNLIKELY(node->depth >= max_depth)) {
        return false;
    }
    
    spatial_hyperoctree_item *items = node->items;
    node->items = NULL;
    node->is_leaf = false;
    
    spatial_num_t cmin[SPATIAL_HYPEROCTREE_DIMS];
    spatial_num_t cmax[SPATIAL_HYPEROCTREE_DIMS];
    
    for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
        spatial_hyperoctree_child_bounds(node, i, cmin, cmax);
        
        node->children[i] = spatial_hyperoctree_node_new(cmin, cmax, node->depth + 1, alloc);
        if (SPATIAL_UNLIKELY(!node->children[i])) {
            for (int j = 0; j < i; j++) {
                spatial_hyperoctree_node_free(node->children[j], alloc, NULL, NULL);
                node->children[j] = NULL;
            }
            node->is_leaf = true;
            node->items = items;
            return false;
        }
    }
    
    /* Redistribute items */
    spatial_hyperoctree_item *item = items;
    while (item) {
        spatial_hyperoctree_item *next = item->next;
        int oct = spatial_hyperoctree_octant(node, item->min, item->max);
        spatial_hyperoctree_node *child = node->children[oct];
        
        item->next = child->items;
        child->items = item;
        child->item_count++;
        
        item = next;
    }
    
    return true;
}

/* ============================================================================
 * Public API
 * ============================================================================ */

SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_new(void) {
    return spatial_hyperoctree_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();
    
    spatial_hyperoctree *ht = (spatial_hyperoctree*)alloc->malloc(sizeof(spatial_hyperoctree), alloc->udata);
    if (SPATIAL_UNLIKELY(!ht)) return NULL;
    
    memset(ht, 0, sizeof(spatial_hyperoctree));
    spatial_refcount_init(&ht->refcount);
    ht->alloc = *alloc;
    ht->max_depth = SPATIAL_HYPEROCTREE_MAX_DEPTH;
    ht->max_items = SPATIAL_HYPEROCTREE_MAX_ITEMS;
    ht->relaxed_atomics = true;
    ht->root = NULL;
    
    return ht;
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_new_with_memento_arena(memento_arena_t *arena) {
    spatial_allocator alloc = spatial_memento_arena_allocator(arena);
    return spatial_hyperoctree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_new_with_memento_pool(memento_pool_t *pool) {
    spatial_allocator alloc = spatial_memento_pool_allocator(pool);
    return spatial_hyperoctree_new_with_allocator(&alloc);
}
#endif

SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_new_with_arena(spatial_arena_t *arena) {
    spatial_allocator alloc = spatial_arena_allocator(arena);
    return spatial_hyperoctree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_new_with_pool(spatial_pool_t *pool) {
    spatial_allocator alloc = spatial_pool_allocator(pool);
    return spatial_hyperoctree_new_with_allocator(&alloc);
}

SPATIAL_INLINE void spatial_hyperoctree_free(spatial_hyperoctree *ht) {
    if (SPATIAL_UNLIKELY(!ht)) return;
    
    if (spatial_refcount_decrement(&ht->refcount) > 0) return;
    
    spatial_hyperoctree_node_free(ht->root, &ht->alloc, &ht->callbacks, ht->udata);
    ht->alloc.free(ht, ht->alloc.udata);
}

SPATIAL_INLINE spatial_hyperoctree* spatial_hyperoctree_clone(const spatial_hyperoctree *ht) {
    if (SPATIAL_UNLIKELY(!ht)) return NULL;
    
    spatial_refcount_increment(&ht->refcount);
    union { const spatial_hyperoctree *in; spatial_hyperoctree *out; } pun; pun.in = ht; return pun.out;
}

SPATIAL_INLINE void spatial_hyperoctree_set_item_callbacks(spatial_hyperoctree *ht,
                                                            const spatial_item_callbacks *cb)
{
    if (ht) ht->callbacks = *cb;
}

SPATIAL_INLINE void spatial_hyperoctree_set_udata(spatial_hyperoctree *ht, void *udata) {
    if (ht) ht->udata = udata;
}

SPATIAL_INLINE void spatial_hyperoctree_opt_relaxed_atomics(spatial_hyperoctree *ht, bool enable) {
    if (ht) ht->relaxed_atomics = enable;
}

SPATIAL_INLINE bool spatial_hyperoctree_init_root(spatial_hyperoctree *ht,
                                                   const spatial_num_t *min,
                                                   const spatial_num_t *max)
{
    if (ht->root) return true;
    
    ht->root = spatial_hyperoctree_node_new(min, max, 0, &ht->alloc);
    return ht->root != NULL;
}


SPATIAL_INLINE bool spatial_hyperoctree_expand_root(spatial_hyperoctree *ht,
                                                     const spatial_num_t *min,
                                                     const spatial_num_t *max)
{
    if (!ht->root) return spatial_hyperoctree_init_root(ht, min, max);
    
    spatial_num_t new_min[SPATIAL_HYPEROCTREE_DIMS];
    spatial_num_t new_max[SPATIAL_HYPEROCTREE_DIMS];
    
    for (int i = 0; i < SPATIAL_HYPEROCTREE_DIMS; i++) {
        new_min[i] = spatial_min(ht->root->min[i], min[i]);
        new_max[i] = spatial_max(ht->root->max[i], max[i]);
        spatial_num_t size = new_max[i] - new_min[i];
        if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
        new_min[i] -= size * (spatial_num_t)0.05;
        new_max[i] += size * (spatial_num_t)0.05;
    }
    
    if (ht->root->is_leaf) {
        for (int i = 0; i < SPATIAL_HYPEROCTREE_DIMS; i++) {
            ht->root->min[i] = new_min[i];
            ht->root->max[i] = new_max[i];
            ht->root->center[i] = (new_min[i] + new_max[i]) * (spatial_num_t)0.5;
        }
        return true;
    }
    
    /* Root has children: create new parent root */
    spatial_hyperoctree_node *new_root = spatial_hyperoctree_node_new(new_min, new_max, 0, &ht->alloc);
    if (!new_root) return false;
    
    new_root->is_leaf = false;
    for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
        new_root->children[i] = NULL;
    }
    
    int old_octant = spatial_hyperoctree_octant(new_root, ht->root->min, ht->root->max);
    new_root->children[old_octant] = ht->root;
    ht->root = new_root;
    return true;
}

SPATIAL_INLINE bool spatial_hyperoctree_insert(spatial_hyperoctree *ht,
                                                const spatial_num_t *SPATIAL_RESTRICT min,
                                                const spatial_num_t *SPATIAL_RESTRICT max,
                                                spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!ht)) return false;
    
    if (SPATIAL_UNLIKELY(!ht->root)) {
        spatial_num_t root_min[SPATIAL_HYPEROCTREE_DIMS];
        spatial_num_t root_max[SPATIAL_HYPEROCTREE_DIMS];
        
        for (int i = 0; i < SPATIAL_HYPEROCTREE_DIMS; i++) {
            root_min[i] = min[i];
            root_max[i] = max[i];
        }
        
        for (int i = 0; i < SPATIAL_HYPEROCTREE_DIMS; i++) {
            spatial_num_t size = root_max[i] - root_min[i];
            if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
            root_min[i] -= size * (spatial_num_t)0.1;
            root_max[i] += size * (spatial_num_t)0.1;
        }
        
        bool ok = spatial_hyperoctree_init_root(ht, root_min, root_max);
        if (!ok) return false;
    }
    
    if (!spatial_bbox_contains(ht->root->min, ht->root->max, min, max, SPATIAL_HYPEROCTREE_DIMS)) {
        if (!spatial_hyperoctree_expand_root(ht, min, max)) return false;
    }
    
    spatial_hyperoctree_item *item = (spatial_hyperoctree_item*)ht->alloc.malloc(
        sizeof(spatial_hyperoctree_item), ht->alloc.udata);
    if (SPATIAL_UNLIKELY(!item)) return false;
    
    spatial_bbox_copy(min, max, item->min, item->max, SPATIAL_HYPEROCTREE_DIMS);
    item->data = data;
    item->next = NULL;
    
    spatial_hyperoctree_node *node = ht->root;
    
    while (!node->is_leaf) {
        int octant = spatial_hyperoctree_octant(node, min, max);
        if (!node->children[octant]) {
            spatial_num_t child_min[SPATIAL_HYPEROCTREE_DIMS];
            spatial_num_t child_max[SPATIAL_HYPEROCTREE_DIMS];
            spatial_hyperoctree_child_bounds(node, octant, child_min, child_max);
            node->children[octant] = spatial_hyperoctree_node_new(child_min, child_max, node->depth + 1, &ht->alloc);
            if (SPATIAL_UNLIKELY(!node->children[octant])) return false;
        }
        node = node->children[octant];
    }
    
    item->next = node->items;
    node->items = item;
    node->item_count++;
    ht->count++;
    
    if (node->item_count > ht->max_items && node->depth < ht->max_depth) {
        spatial_hyperoctree_node_split(node, ht->max_depth, &ht->alloc);
    }
    
    return true;
}

SPATIAL_INLINE void spatial_hyperoctree_search(const spatial_hyperoctree *ht,
                                                const spatial_num_t *SPATIAL_RESTRICT qmin,
                                                const spatial_num_t *SPATIAL_RESTRICT qmax,
                                                spatial_hyperoctree_iter_fn iter,
                                                void *udata)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root || !iter)) return;

    spatial_hyperoctree_node *stk[(SPATIAL_HYPEROCTREE_MAX_DEPTH + 1) * SPATIAL_HYPEROCTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = ht->root;

    while (sp > 0) {
        spatial_hyperoctree_node *node = stk[--sp];

        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(node->min, node->max, qmin, qmax, SPATIAL_HYPEROCTREE_DIMS))) {
            continue;
        }

        if (node->is_leaf) {
            spatial_hyperoctree_item *item = node->items;
            while (item) {
                if (spatial_bbox_overlaps(item->min, item->max, qmin, qmax, SPATIAL_HYPEROCTREE_DIMS)) {
                    if (!iter(item->min, item->max, item->data, udata)) {
                        return;
                    }
                }
                item = item->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE void spatial_hyperoctree_scan(const spatial_hyperoctree *ht,
                                              spatial_hyperoctree_iter_fn iter,
                                              void *udata)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root || !iter)) return;

    spatial_hyperoctree_node *stk[(SPATIAL_HYPEROCTREE_MAX_DEPTH + 1) * SPATIAL_HYPEROCTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = ht->root;

    while (sp > 0) {
        spatial_hyperoctree_node *node = stk[--sp];

        if (node->is_leaf) {
            spatial_hyperoctree_item *item = node->items;
            while (item) {
                if (!iter(item->min, item->max, item->data, udata)) {
                    return;
                }
                item = item->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE int spatial_hyperoctree_count(const spatial_hyperoctree *ht) {
    return ht ? ht->count : 0;
}

SPATIAL_INLINE bool spatial_hyperoctree_delete(spatial_hyperoctree *ht,
                                                const spatial_num_t *min,
                                                const spatial_num_t *max,
                                                spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root)) return false;

    spatial_hyperoctree_node *stk[(SPATIAL_HYPEROCTREE_MAX_DEPTH + 1) * SPATIAL_HYPEROCTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = ht->root;

    while (sp > 0 && !found) {
        spatial_hyperoctree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_HYPEROCTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_hyperoctree_item **pp = &node->items;
            while (*pp) {
                spatial_hyperoctree_item *item = *pp;
                if (item->data == data &&
                    spatial_bbox_overlaps(item->min, item->max, min, max, SPATIAL_HYPEROCTREE_DIMS)) {
                    *pp = item->next;
                    if (ht->callbacks.free) {
                        ht->callbacks.free(item->data, ht->udata);
                    }
                    ht->alloc.free(item, ht->alloc.udata);
                    node->item_count--;
                    ht->count--;
                    found = true;
                    break;
                }
                pp = &(*pp)->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }

    return found;
}

typedef bool (*spatial_hyperoctree_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_hyperoctree_delete_with_comparator(spatial_hyperoctree *ht,
                                                                const spatial_num_t *min,
                                                                const spatial_num_t *max,
                                                                spatial_hyperoctree_cmp_fn cmp,
                                                                void *cmp_udata)
{
    if (SPATIAL_UNLIKELY(!ht || !ht->root)) return false;

    spatial_hyperoctree_node *stk[(SPATIAL_HYPEROCTREE_MAX_DEPTH + 1) * SPATIAL_HYPEROCTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = ht->root;

    while (sp > 0 && !found) {
        spatial_hyperoctree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_HYPEROCTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_hyperoctree_item **pp = &node->items;
            while (*pp) {
                spatial_hyperoctree_item *item = *pp;
                if (cmp(item->data, cmp_udata) &&
                    spatial_bbox_overlaps(item->min, item->max, min, max, SPATIAL_HYPEROCTREE_DIMS)) {
                    *pp = item->next;
                    if (ht->callbacks.free) {
                        ht->callbacks.free(item->data, ht->udata);
                    }
                    ht->alloc.free(item, ht->alloc.udata);
                    node->item_count--;
                    ht->count--;
                    found = true;
                    break;
                }
                pp = &(*pp)->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_HYPEROCTREE_CHILDREN; i++) {
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

#endif /* SPATIAL_HYPEROCTREE_H */
