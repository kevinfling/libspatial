/**
 * @file   octree.h
 * @brief  3D Octree for spatial indexing
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * An octree recursively subdivides 3D space into eight octants.
 * Used for 3D collision detection, point cloud indexing, and volumetric queries.
 */

#ifndef SPATIAL_OCTREE_H
#define SPATIAL_OCTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

#define SPATIAL_OCTREE_DIMS 3
#define SPATIAL_OCTREE_CHILDREN 8

#ifndef SPATIAL_OCTREE_MAX_DEPTH
    #define SPATIAL_OCTREE_MAX_DEPTH 16
#endif

#ifndef SPATIAL_OCTREE_MAX_ITEMS
    #define SPATIAL_OCTREE_MAX_ITEMS 16
#endif

/* Forward declarations */
struct spatial_octree_node;
struct spatial_octree_item;

/* Item structure */
typedef struct spatial_octree_item {
    spatial_num_t min[SPATIAL_OCTREE_DIMS];
    spatial_num_t max[SPATIAL_OCTREE_DIMS];
    spatial_data_t data;
    struct spatial_octree_item *next;
} spatial_octree_item;

/* Node structure (cache-line aligned for SIMD-ready traversal) */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_octree_node {
    spatial_num_t min[SPATIAL_OCTREE_DIMS];
    spatial_num_t max[SPATIAL_OCTREE_DIMS];
    spatial_num_t center[SPATIAL_OCTREE_DIMS];
    
    union {
        struct spatial_octree_node *children[SPATIAL_OCTREE_CHILDREN];
        spatial_octree_item *items;
    };
    
    int item_count;
    int depth;
    bool is_leaf;
} spatial_octree_node;

/* Main octree structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_octree {
    spatial_octree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    int count;
    int max_depth;
    int max_items;
    bool relaxed_atomics;
} spatial_octree;

/* Iterator callback */
typedef bool (*spatial_octree_iter_fn)(const spatial_num_t *min,
                                        const spatial_num_t *max,
                                        spatial_data_t data,
                                        void *udata);

/* Forward declarations for public API */
SPATIAL_INLINE spatial_octree* spatial_octree_new_with_allocator(const spatial_allocator *alloc);

/* ============================================================================
 * Internal Node Management
 * ============================================================================ */

SPATIAL_INLINE spatial_octree_node* spatial_octree_node_new(
    const spatial_num_t *min,
    const spatial_num_t *max,
    int depth,
    const spatial_allocator *alloc)
{
    spatial_octree_node *node = (spatial_octree_node*)alloc->malloc(sizeof(spatial_octree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    spatial_bbox_copy(min, max, node->min, node->max, SPATIAL_OCTREE_DIMS);
    
    for (int i = 0; i < SPATIAL_OCTREE_DIMS; i++) {
        node->center[i] = (min[i] + max[i]) * (spatial_num_t)0.5;
    }
    
    for (int i = 0; i < SPATIAL_OCTREE_CHILDREN; i++) {
        node->children[i] = NULL;
    }
    
    node->item_count = 0;
    node->depth = depth;
    node->is_leaf = true;
    
    return node;
}

SPATIAL_INLINE void spatial_octree_node_free(spatial_octree_node *node,
                                              const spatial_allocator *alloc,
                                              const spatial_item_callbacks *cb,
                                              void *udata);

SPATIAL_INLINE void spatial_octree_item_free_list(spatial_octree_item *item,
                                                   const spatial_allocator *alloc,
                                                   const spatial_item_callbacks *cb,
                                                   void *udata)
{
    while (item) {
        spatial_octree_item *next = item->next;
        if (cb && cb->free) {
            cb->free(item->data, udata);
        }
        alloc->free(item, alloc->udata);
        item = next;
    }
}

SPATIAL_INLINE void spatial_octree_node_free(spatial_octree_node *node,
                                              const spatial_allocator *alloc,
                                              const spatial_item_callbacks *cb,
                                              void *udata)
{
    if (SPATIAL_UNLIKELY(!node)) return;
    
    spatial_octree_node *stk[SPATIAL_OCTREE_MAX_DEPTH * SPATIAL_OCTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = node;
    
    while (sp > 0) {
        spatial_octree_node *n = stk[--sp];
        
        if (n->is_leaf) {
            spatial_octree_item_free_list(n->items, alloc, cb, udata);
        } else {
            for (int i = 0; i < SPATIAL_OCTREE_CHILDREN; i++) {
                if (n->children[i]) stk[sp++] = n->children[i];
            }
        }
        alloc->free(n, alloc->udata);
    }
}

SPATIAL_INLINE int spatial_octree_octant(const spatial_octree_node *node,
                                          const spatial_num_t *min,
                                          const spatial_num_t *max)
{
    spatial_num_t px = (min[0] + max[0]) * (spatial_num_t)0.5;
    spatial_num_t py = (min[1] + max[1]) * (spatial_num_t)0.5;
    spatial_num_t pz = (min[2] + max[2]) * (spatial_num_t)0.5;
    
    return ((px >= node->center[0]) ? 1 : 0)
         | ((py >= node->center[1]) ? 2 : 0)
         | ((pz >= node->center[2]) ? 4 : 0);
}

SPATIAL_INLINE void spatial_octree_child_bounds(const spatial_octree_node *node,
                                                 int octant,
                                                 spatial_num_t *child_min,
                                                 spatial_num_t *child_max)
{
    child_min[0] = (octant & 1) ? node->center[0] : node->min[0];
    child_max[0] = (octant & 1) ? node->max[0] : node->center[0];
    child_min[1] = (octant & 2) ? node->center[1] : node->min[1];
    child_max[1] = (octant & 2) ? node->max[1] : node->center[1];
    child_min[2] = (octant & 4) ? node->center[2] : node->min[2];
    child_max[2] = (octant & 4) ? node->max[2] : node->center[2];
}

SPATIAL_INLINE bool spatial_octree_node_split(spatial_octree_node *node,
                                               const spatial_allocator *alloc)
{
    if (SPATIAL_UNLIKELY(node->depth >= SPATIAL_OCTREE_MAX_DEPTH)) {
        return false;
    }
    
    spatial_octree_item *items = node->items;
    node->items = NULL;
    node->is_leaf = false;
    
    /* Redistribute items, allocating children on demand */
    spatial_octree_item *item = items;
    while (item) {
        spatial_octree_item *next = item->next;
        int octant = spatial_octree_octant(node, item->min, item->max);
        
        if (SPATIAL_UNLIKELY(!node->children[octant])) {
            spatial_num_t cmin[SPATIAL_OCTREE_DIMS], cmax[SPATIAL_OCTREE_DIMS];
            spatial_octree_child_bounds(node, octant, cmin, cmax);
            node->children[octant] = spatial_octree_node_new(cmin, cmax, node->depth + 1, alloc);
            if (SPATIAL_UNLIKELY(!node->children[octant])) {
                for (int j = 0; j < SPATIAL_OCTREE_CHILDREN; j++) {
                    if (node->children[j]) {
                        alloc->free(node->children[j], alloc->udata);
                        node->children[j] = NULL;
                    }
                }
                node->is_leaf = true;
                node->items = items;
                return false;
            }
        }
        
        spatial_octree_node *child = node->children[octant];
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

SPATIAL_INLINE spatial_octree* spatial_octree_new(void) {
    return spatial_octree_new_with_allocator(NULL);
}

SPATIAL_INLINE spatial_octree* spatial_octree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();
    
    spatial_octree *ot = (spatial_octree*)alloc->malloc(sizeof(spatial_octree), alloc->udata);
    if (SPATIAL_UNLIKELY(!ot)) return NULL;
    
    memset(ot, 0, sizeof(spatial_octree));
    spatial_refcount_init(&ot->refcount);
    ot->alloc = *alloc;
    ot->max_depth = SPATIAL_OCTREE_MAX_DEPTH;
    ot->max_items = SPATIAL_OCTREE_MAX_ITEMS;
    ot->relaxed_atomics = true;
    ot->root = NULL;
    
    return ot;
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_octree* spatial_octree_new_with_memento_arena(memento_arena_t *arena) {
    spatial_allocator alloc = spatial_memento_arena_allocator(arena);
    return spatial_octree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_octree* spatial_octree_new_with_memento_pool(memento_pool_t *pool) {
    spatial_allocator alloc = spatial_memento_pool_allocator(pool);
    return spatial_octree_new_with_allocator(&alloc);
}
#endif

SPATIAL_INLINE spatial_octree* spatial_octree_new_with_arena(spatial_arena_t *arena) {
    spatial_allocator alloc = spatial_arena_allocator(arena);
    return spatial_octree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_octree* spatial_octree_new_with_pool(spatial_pool_t *pool) {
    spatial_allocator alloc = spatial_pool_allocator(pool);
    return spatial_octree_new_with_allocator(&alloc);
}



SPATIAL_INLINE void spatial_octree_free(spatial_octree *ot) {
    if (SPATIAL_UNLIKELY(!ot)) return;
    
    if (spatial_refcount_decrement(&ot->refcount) > 0) return;
    
    spatial_octree_node_free(ot->root, &ot->alloc, &ot->callbacks, ot->udata);
    ot->alloc.free(ot, ot->alloc.udata);
}

SPATIAL_INLINE spatial_octree* spatial_octree_clone(const spatial_octree *ot) {
    if (SPATIAL_UNLIKELY(!ot)) return NULL;
    
    spatial_refcount_increment(&ot->refcount);
    union { const spatial_octree *in; spatial_octree *out; } pun; pun.in = ot; return pun.out;
}

SPATIAL_INLINE void spatial_octree_set_item_callbacks(spatial_octree *ot,
                                                       const spatial_item_callbacks *cb)
{
    if (ot) ot->callbacks = *cb;
}

SPATIAL_INLINE void spatial_octree_set_udata(spatial_octree *ot, void *udata) {
    if (ot) ot->udata = udata;
}

SPATIAL_INLINE void spatial_octree_opt_relaxed_atomics(spatial_octree *ot, bool enable) {
    if (ot) ot->relaxed_atomics = enable;
}

SPATIAL_INLINE bool spatial_octree_init_root(spatial_octree *ot,
                                              const spatial_num_t *min,
                                              const spatial_num_t *max)
{
    if (ot->root) return true;
    
    ot->root = spatial_octree_node_new(min, max, 0, &ot->alloc);
    return ot->root != NULL;
}


SPATIAL_INLINE bool spatial_octree_expand_root(spatial_octree *ot,
                                                const spatial_num_t *min,
                                                const spatial_num_t *max)
{
    if (!ot->root) return spatial_octree_init_root(ot, min, max);
    
    spatial_num_t new_min[SPATIAL_OCTREE_DIMS];
    spatial_num_t new_max[SPATIAL_OCTREE_DIMS];
    
    for (int i = 0; i < SPATIAL_OCTREE_DIMS; i++) {
        new_min[i] = spatial_min(ot->root->min[i], min[i]);
        new_max[i] = spatial_max(ot->root->max[i], max[i]);
        spatial_num_t size = new_max[i] - new_min[i];
        if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
        new_min[i] -= size * (spatial_num_t)0.05;
        new_max[i] += size * (spatial_num_t)0.05;
    }
    
    if (ot->root->is_leaf) {
        for (int i = 0; i < SPATIAL_OCTREE_DIMS; i++) {
            ot->root->min[i] = new_min[i];
            ot->root->max[i] = new_max[i];
            ot->root->center[i] = (new_min[i] + new_max[i]) * (spatial_num_t)0.5;
        }
        return true;
    }
    
    /* Root has children: create new parent root */
    spatial_octree_node *new_root = spatial_octree_node_new(new_min, new_max, 0, &ot->alloc);
    if (!new_root) return false;
    
    new_root->is_leaf = false;
    for (int i = 0; i < 8; i++) {
        new_root->children[i] = NULL;
    }
    
    int old_octant = spatial_octree_octant(new_root, ot->root->min, ot->root->max);
    new_root->children[old_octant] = ot->root;
    ot->root = new_root;
    return true;
}

SPATIAL_INLINE bool spatial_octree_insert(spatial_octree *ot,
                                           const spatial_num_t *SPATIAL_RESTRICT min,
                                           const spatial_num_t *SPATIAL_RESTRICT max,
                                           spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!ot)) return false;
    
    if (SPATIAL_UNLIKELY(!ot->root)) {
        spatial_num_t root_min[SPATIAL_OCTREE_DIMS];
        spatial_num_t root_max[SPATIAL_OCTREE_DIMS];
        
        for (int i = 0; i < SPATIAL_OCTREE_DIMS; i++) {
            root_min[i] = min[i];
            root_max[i] = max[i];
        }
        
        for (int i = 0; i < SPATIAL_OCTREE_DIMS; i++) {
            spatial_num_t size = root_max[i] - root_min[i];
            if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
            root_min[i] -= size * (spatial_num_t)0.1;
            root_max[i] += size * (spatial_num_t)0.1;
        }
        
        if (!spatial_octree_init_root(ot, root_min, root_max)) return false;
    }
    
    if (!spatial_bbox_contains(ot->root->min, ot->root->max, min, max, SPATIAL_OCTREE_DIMS)) {
        if (!spatial_octree_expand_root(ot, min, max)) return false;
    }
    
    spatial_octree_item *item = (spatial_octree_item*)ot->alloc.malloc(sizeof(spatial_octree_item), ot->alloc.udata);
    if (SPATIAL_UNLIKELY(!item)) return false;
    
    spatial_bbox_copy(min, max, item->min, item->max, SPATIAL_OCTREE_DIMS);
    item->data = data;
    item->next = NULL;
    
    spatial_octree_node *node = ot->root;
    
    while (!node->is_leaf) {
        int octant = spatial_octree_octant(node, min, max);
        if (!node->children[octant]) {
            spatial_num_t child_min[SPATIAL_OCTREE_DIMS];
            spatial_num_t child_max[SPATIAL_OCTREE_DIMS];
            spatial_octree_child_bounds(node, octant, child_min, child_max);
            node->children[octant] = spatial_octree_node_new(child_min, child_max, node->depth + 1, &ot->alloc);
            if (SPATIAL_UNLIKELY(!node->children[octant])) return false;
        }
        node = node->children[octant];
    }
    
    item->next = node->items;
    node->items = item;
    node->item_count++;
    ot->count++;
    
    if (node->item_count > ot->max_items && node->depth < ot->max_depth) {
        spatial_octree_node_split(node, &ot->alloc);
    }
    
    return true;
}

SPATIAL_INLINE void spatial_octree_search(const spatial_octree *ot,
                                           const spatial_num_t *SPATIAL_RESTRICT qmin,
                                           const spatial_num_t *SPATIAL_RESTRICT qmax,
                                           spatial_octree_iter_fn iter,
                                           void *udata)
{
    if (SPATIAL_UNLIKELY(!ot || !ot->root || !iter)) return;

    spatial_octree_node *stk[(SPATIAL_OCTREE_MAX_DEPTH + 1) * SPATIAL_OCTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = ot->root;

    while (sp > 0) {
        spatial_octree_node *node = stk[--sp];

        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(node->min, node->max, qmin, qmax, SPATIAL_OCTREE_DIMS))) {
            continue;
        }

        if (node->is_leaf) {
            spatial_octree_item *item = node->items;
            while (item) {
                if (spatial_bbox_overlaps(item->min, item->max, qmin, qmax, SPATIAL_OCTREE_DIMS)) {
                    if (!iter(item->min, item->max, item->data, udata)) {
                        return;
                    }
                }
                item = item->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_OCTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    SPATIAL_PREFETCH(node->children[i], 0, 3);
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE void spatial_octree_scan(const spatial_octree *ot,
                                         spatial_octree_iter_fn iter,
                                         void *udata)
{
    if (SPATIAL_UNLIKELY(!ot || !ot->root || !iter)) return;

    spatial_octree_node *stk[(SPATIAL_OCTREE_MAX_DEPTH + 1) * SPATIAL_OCTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = ot->root;

    while (sp > 0) {
        spatial_octree_node *node = stk[--sp];

        if (node->is_leaf) {
            spatial_octree_item *item = node->items;
            while (item) {
                if (!iter(item->min, item->max, item->data, udata)) {
                    return;
                }
                item = item->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_OCTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

SPATIAL_INLINE int spatial_octree_count(const spatial_octree *ot) {
    return ot ? ot->count : 0;
}

SPATIAL_INLINE bool spatial_octree_delete(spatial_octree *ot,
                                           const spatial_num_t *min,
                                           const spatial_num_t *max,
                                           spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!ot || !ot->root)) return false;

    spatial_octree_node *stk[(SPATIAL_OCTREE_MAX_DEPTH + 1) * SPATIAL_OCTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = ot->root;

    while (sp > 0 && !found) {
        spatial_octree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_OCTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_octree_item **pp = &node->items;
            while (*pp) {
                spatial_octree_item *item = *pp;
                if (item->data == data &&
                    spatial_bbox_overlaps(item->min, item->max, min, max, SPATIAL_OCTREE_DIMS)) {
                    *pp = item->next;
                    if (ot->callbacks.free) {
                        ot->callbacks.free(item->data, ot->udata);
                    }
                    ot->alloc.free(item, ot->alloc.udata);
                    node->item_count--;
                    ot->count--;
                    found = true;
                    break;
                }
                pp = &(*pp)->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_OCTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    stk[sp++] = node->children[i];
                }
            }
        }
    }

    return found;
}

typedef bool (*spatial_octree_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_octree_delete_with_comparator(spatial_octree *ot,
                                                           const spatial_num_t *min,
                                                           const spatial_num_t *max,
                                                           spatial_octree_cmp_fn cmp,
                                                           void *cmp_udata)
{
    if (SPATIAL_UNLIKELY(!ot || !ot->root)) return false;

    spatial_octree_node *stk[(SPATIAL_OCTREE_MAX_DEPTH + 1) * SPATIAL_OCTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = ot->root;

    while (sp > 0 && !found) {
        spatial_octree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_OCTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_octree_item **pp = &node->items;
            while (*pp) {
                spatial_octree_item *item = *pp;
                if (cmp(item->data, cmp_udata) &&
                    spatial_bbox_overlaps(item->min, item->max, min, max, SPATIAL_OCTREE_DIMS)) {
                    *pp = item->next;
                    if (ot->callbacks.free) {
                        ot->callbacks.free(item->data, ot->udata);
                    }
                    ot->alloc.free(item, ot->alloc.udata);
                    node->item_count--;
                    ot->count--;
                    found = true;
                    break;
                }
                pp = &(*pp)->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_OCTREE_CHILDREN; i++) {
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

#endif /* SPATIAL_OCTREE_H */
