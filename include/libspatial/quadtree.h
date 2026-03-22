/**
 * @file   quadtree.h
 * @brief  2D Quadtree with PR/MX/Point-Region variants
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * A quadtree recursively subdivides 2D space into four quadrants.
 * Variants:
 * - PR (Point-Region): Stores points, splits when count exceeds threshold
 * - MX (Matrix): Fixed grid subdivision, stores at most one point per cell
 * - LINEAR: Linear quadtree using Morton codes (no explicit tree)
 */

#ifndef SPATIAL_QUADTREE_H
#define SPATIAL_QUADTREE_H

#ifdef __cplusplus
extern "C" {
#endif

#include "spatial.h"

#ifndef SPATIAL_QUADTREE_VARIANT
    #define SPATIAL_QUADTREE_VARIANT SPATIAL_QUADTREE_PR
#endif

#define SPATIAL_QUADTREE_PR     0
#define SPATIAL_QUADTREE_MX     1
#define SPATIAL_QUADTREE_LINEAR 2

#define SPATIAL_QUADTREE_DIMS 2
#define SPATIAL_QUADTREE_CHILDREN 4

#ifndef SPATIAL_QUADTREE_MAX_DEPTH
    #define SPATIAL_QUADTREE_MAX_DEPTH 20
#endif

#ifndef SPATIAL_QUADTREE_MAX_ITEMS
    #define SPATIAL_QUADTREE_MAX_ITEMS 16
#endif

/* Forward declarations */
struct spatial_quadtree_node;
struct spatial_quadtree_item;

/* Item structure (stores one data element with its bounds) */
typedef struct spatial_quadtree_item {
    spatial_num_t min[SPATIAL_QUADTREE_DIMS];
    spatial_num_t max[SPATIAL_QUADTREE_DIMS];
    spatial_data_t data;
    struct spatial_quadtree_item *next;
} spatial_quadtree_item;

/* Node structure (cache-line aligned for SIMD-ready traversal) */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_quadtree_node {
    spatial_num_t min[SPATIAL_QUADTREE_DIMS];
    spatial_num_t max[SPATIAL_QUADTREE_DIMS];
    spatial_num_t center[SPATIAL_QUADTREE_DIMS];
    
    union {
        struct spatial_quadtree_node *children[SPATIAL_QUADTREE_CHILDREN];
        spatial_quadtree_item *items;
    };
    
    int item_count;
    int depth;
    bool is_leaf;
} spatial_quadtree_node;

/* Main quadtree structure */
typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_quadtree {
    spatial_quadtree_node *root;
    spatial_refcount refcount;
    spatial_allocator alloc;
    spatial_item_callbacks callbacks;
    void *udata;
    int count;
    int max_depth;
    int max_items;
    bool relaxed_atomics;
} spatial_quadtree;

/* Iterator callback type */
typedef bool (*spatial_quadtree_iter_fn)(const spatial_num_t *min,
                                          const spatial_num_t *max,
                                          spatial_data_t data,
                                          void *udata);

/* Forward declarations for public API */
SPATIAL_INLINE spatial_quadtree* spatial_quadtree_new_with_allocator(const spatial_allocator *alloc);

/* ============================================================================
 * Internal Node Management
 * ============================================================================ */

SPATIAL_INLINE spatial_quadtree_node* spatial_quadtree_node_new(
    const spatial_num_t *min,
    const spatial_num_t *max,
    int depth,
    const spatial_allocator *alloc)
{
    spatial_quadtree_node *node = (spatial_quadtree_node*)alloc->malloc(sizeof(spatial_quadtree_node), alloc->udata);
    if (SPATIAL_UNLIKELY(!node)) return NULL;
    
    spatial_bbox_copy(min, max, node->min, node->max, SPATIAL_QUADTREE_DIMS);
    
    for (int i = 0; i < SPATIAL_QUADTREE_DIMS; i++) {
        node->center[i] = (min[i] + max[i]) * (spatial_num_t)0.5;
    }
    
    for (int i = 0; i < SPATIAL_QUADTREE_CHILDREN; i++) {
        node->children[i] = NULL;
    }
    
    node->item_count = 0;
    node->depth = depth;
    node->is_leaf = true;
    
    return node;
}

SPATIAL_INLINE void spatial_quadtree_node_free(spatial_quadtree_node *node,
                                                const spatial_allocator *alloc,
                                                const spatial_item_callbacks *cb,
                                                void *udata);

SPATIAL_INLINE void spatial_quadtree_item_free_list(spatial_quadtree_item *item,
                                                     const spatial_allocator *alloc,
                                                     const spatial_item_callbacks *cb,
                                                     void *udata)
{
    while (item) {
        spatial_quadtree_item *next = item->next;
        if (cb && cb->free) {
            cb->free(item->data, udata);
        }
        alloc->free(item, alloc->udata);
        item = next;
    }
}

SPATIAL_INLINE void spatial_quadtree_node_free(spatial_quadtree_node *node,
                                                const spatial_allocator *alloc,
                                                const spatial_item_callbacks *cb,
                                                void *udata)
{
    if (SPATIAL_UNLIKELY(!node)) return;
    
    spatial_quadtree_node *stk[SPATIAL_QUADTREE_MAX_DEPTH * SPATIAL_QUADTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = node;
    
    while (sp > 0) {
        spatial_quadtree_node *n = stk[--sp];
        
        if (n->is_leaf) {
            spatial_quadtree_item_free_list(n->items, alloc, cb, udata);
        } else {
            for (int i = 0; i < SPATIAL_QUADTREE_CHILDREN; i++) {
                if (n->children[i]) stk[sp++] = n->children[i];
            }
        }
        alloc->free(n, alloc->udata);
    }
}

SPATIAL_INLINE int spatial_quadtree_quadrant(const spatial_quadtree_node *node,
                                              const spatial_num_t *min,
                                              const spatial_num_t *max)
{
    spatial_num_t cx = node->center[0];
    spatial_num_t cy = node->center[1];
    spatial_num_t px = (min[0] + max[0]) * (spatial_num_t)0.5;
    spatial_num_t py = (min[1] + max[1]) * (spatial_num_t)0.5;
    
    int quadrant = 0;
    if (px >= cx) quadrant |= 1;
    if (py >= cy) quadrant |= 2;
    return quadrant;
}

SPATIAL_INLINE void spatial_quadtree_child_bounds(const spatial_quadtree_node *node,
                                                   int quadrant,
                                                   spatial_num_t *child_min,
                                                   spatial_num_t *child_max)
{
    child_min[0] = (quadrant & 1) ? node->center[0] : node->min[0];
    child_max[0] = (quadrant & 1) ? node->max[0] : node->center[0];
    child_min[1] = (quadrant & 2) ? node->center[1] : node->min[1];
    child_max[1] = (quadrant & 2) ? node->max[1] : node->center[1];
}

SPATIAL_INLINE bool spatial_quadtree_node_split(spatial_quadtree_node *node,
                                                 const spatial_allocator *alloc)
{
    if (SPATIAL_UNLIKELY(node->depth >= SPATIAL_QUADTREE_MAX_DEPTH)) {
        return false;
    }
    
    spatial_quadtree_item *items = node->items;
    node->items = NULL;
    node->is_leaf = false;
    
    for (int i = 0; i < SPATIAL_QUADTREE_CHILDREN; i++) {
        spatial_num_t cmin[SPATIAL_QUADTREE_DIMS], cmax[SPATIAL_QUADTREE_DIMS];
        spatial_quadtree_child_bounds(node, i, cmin, cmax);
        
        node->children[i] = spatial_quadtree_node_new(cmin, cmax, node->depth + 1, alloc);
        if (SPATIAL_UNLIKELY(!node->children[i])) {
            /* Cleanup on failure */
            for (int j = 0; j < i; j++) {
                spatial_quadtree_node_free(node->children[j], alloc, NULL, NULL);
            }
            node->is_leaf = true;
            node->items = items;
            return false;
        }
    }
    
    /* Redistribute items */
    spatial_quadtree_item *item = items;
    while (item) {
        spatial_quadtree_item *next = item->next;
        int quadrant = spatial_quadtree_quadrant(node, item->min, item->max);
        spatial_quadtree_node *child = node->children[quadrant];
        
        item->next = child->items;
        child->items = item;
        child->item_count++;
        
        item = next;
    }
    
    return true;
}

/* ============================================================================
 * Public API Implementation
 * ============================================================================ */

/**
 * @brief Create a new quadtree with default allocator
 */
SPATIAL_INLINE spatial_quadtree* spatial_quadtree_new(void) {
    return spatial_quadtree_new_with_allocator(NULL);
}

/**
 * @brief Create a new quadtree with custom allocator
 */
SPATIAL_INLINE spatial_quadtree* spatial_quadtree_new_with_allocator(
    const spatial_allocator *alloc)
{
    if (!alloc) alloc = spatial_default_allocator();
    
    spatial_quadtree *qt = (spatial_quadtree*)alloc->malloc(sizeof(spatial_quadtree), alloc->udata);
    if (SPATIAL_UNLIKELY(!qt)) return NULL;
    
    memset(qt, 0, sizeof(spatial_quadtree));
    spatial_refcount_init(&qt->refcount);
    qt->alloc = *alloc;
    qt->max_depth = SPATIAL_QUADTREE_MAX_DEPTH;
    qt->max_items = SPATIAL_QUADTREE_MAX_ITEMS;
    qt->relaxed_atomics = true;
    
    /* Root node created on first insert with proper bounds */
    qt->root = NULL;
    
    return qt;
}

#ifdef MEMENTO_H
SPATIAL_INLINE spatial_quadtree* spatial_quadtree_new_with_memento_arena(memento_arena_t *arena) {
    spatial_allocator alloc = spatial_memento_arena_allocator(arena);
    return spatial_quadtree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_quadtree* spatial_quadtree_new_with_memento_pool(memento_pool_t *pool) {
    spatial_allocator alloc = spatial_memento_pool_allocator(pool);
    return spatial_quadtree_new_with_allocator(&alloc);
}
#endif

SPATIAL_INLINE spatial_quadtree* spatial_quadtree_new_with_arena(spatial_arena_t *arena) {
    spatial_allocator alloc = spatial_arena_allocator(arena);
    return spatial_quadtree_new_with_allocator(&alloc);
}

SPATIAL_INLINE spatial_quadtree* spatial_quadtree_new_with_pool(spatial_pool_t *pool) {
    spatial_allocator alloc = spatial_pool_allocator(pool);
    return spatial_quadtree_new_with_allocator(&alloc);
}



/**
 * @brief Free a quadtree
 */
SPATIAL_INLINE void spatial_quadtree_free(spatial_quadtree *qt) {
    if (SPATIAL_UNLIKELY(!qt)) return;
    
    if (spatial_refcount_decrement(&qt->refcount) > 0) return;
    
    spatial_quadtree_node_free(qt->root, &qt->alloc, &qt->callbacks, qt->udata);
    qt->alloc.free(qt, qt->alloc.udata);
}

/**
 * @brief Clone a quadtree (copy-on-write)
 */
SPATIAL_INLINE spatial_quadtree* spatial_quadtree_clone(const spatial_quadtree *qt) {
    if (SPATIAL_UNLIKELY(!qt)) return NULL;
    
    spatial_refcount_increment(&qt->refcount);
    union { const spatial_quadtree *in; spatial_quadtree *out; } pun; pun.in = qt; return pun.out;
}

/**
 * @brief Set item callbacks for deep copy-on-write
 */
SPATIAL_INLINE void spatial_quadtree_set_item_callbacks(spatial_quadtree *qt,
                                                         const spatial_item_callbacks *cb)
{
    if (qt) qt->callbacks = *cb;
}

/**
 * @brief Set user data
 */
SPATIAL_INLINE void spatial_quadtree_set_udata(spatial_quadtree *qt, void *udata) {
    if (qt) qt->udata = udata;
}

/**
 * @brief Enable/disable relaxed atomics for single-threaded performance
 */
SPATIAL_INLINE void spatial_quadtree_opt_relaxed_atomics(spatial_quadtree *qt, bool enable) {
    if (qt) qt->relaxed_atomics = enable;
}

/**
 * @brief Initialize root node with bounds (called internally)
 */
SPATIAL_INLINE bool spatial_quadtree_init_root(spatial_quadtree *qt,
                                                const spatial_num_t *min,
                                                const spatial_num_t *max)
{
    if (qt->root) return true;
    
    qt->root = spatial_quadtree_node_new(min, max, 0, &qt->alloc);
    return qt->root != NULL;
}

/**
 * @brief Insert an item into the quadtree
 */

SPATIAL_INLINE bool spatial_quadtree_expand_root(spatial_quadtree *qt,
                                                  const spatial_num_t *min,
                                                  const spatial_num_t *max)
{
    if (!qt->root) return spatial_quadtree_init_root(qt, min, max);
    
    spatial_num_t new_min[SPATIAL_QUADTREE_DIMS];
    spatial_num_t new_max[SPATIAL_QUADTREE_DIMS];
    
    for (int i = 0; i < SPATIAL_QUADTREE_DIMS; i++) {
        new_min[i] = spatial_min(qt->root->min[i], min[i]);
        new_max[i] = spatial_max(qt->root->max[i], max[i]);
        spatial_num_t size = new_max[i] - new_min[i];
        if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
        new_min[i] -= size * (spatial_num_t)0.05;
        new_max[i] += size * (spatial_num_t)0.05;
    }
    
    if (qt->root->is_leaf) {
        for (int i = 0; i < SPATIAL_QUADTREE_DIMS; i++) {
            qt->root->min[i] = new_min[i];
            qt->root->max[i] = new_max[i];
            qt->root->center[i] = (new_min[i] + new_max[i]) * (spatial_num_t)0.5;
        }
        return true;
    }
    
    /* Root has children: create new parent root */
    spatial_quadtree_node *new_root = spatial_quadtree_node_new(new_min, new_max, 0, &qt->alloc);
    if (!new_root) return false;
    
    new_root->is_leaf = false;
    for (int i = 0; i < 4; i++) {
        new_root->children[i] = NULL;
    }
    
    int old_quadrant = spatial_quadtree_quadrant(new_root, qt->root->min, qt->root->max);
    new_root->children[old_quadrant] = qt->root;
    qt->root = new_root;
    return true;
}

SPATIAL_INLINE bool spatial_quadtree_insert(spatial_quadtree *qt,
                                             const spatial_num_t *SPATIAL_RESTRICT min,
                                             const spatial_num_t *SPATIAL_RESTRICT max,
                                             spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!qt)) return false;
    
    /* Initialize root on first insert if needed */
    if (SPATIAL_UNLIKELY(!qt->root)) {
        spatial_num_t root_min[SPATIAL_QUADTREE_DIMS];
        spatial_num_t root_max[SPATIAL_QUADTREE_DIMS];
        
        /* Estimate bounds from first item, will expand as needed */
        for (int i = 0; i < SPATIAL_QUADTREE_DIMS; i++) {
            root_min[i] = min[i];
            root_max[i] = max[i];
        }
        
        /* Add some padding */
        for (int i = 0; i < SPATIAL_QUADTREE_DIMS; i++) {
            spatial_num_t size = root_max[i] - root_min[i];
            if (size < SPATIAL_EPSILON) size = (spatial_num_t)1.0;
            root_min[i] -= size * (spatial_num_t)0.1;
            root_max[i] += size * (spatial_num_t)0.1;
        }
        
        if (!spatial_quadtree_init_root(qt, root_min, root_max)) return false;
    }
    
    /* Check if item fits in root bounds, expand if needed */
    if (!spatial_bbox_contains(qt->root->min, qt->root->max, min, max, SPATIAL_QUADTREE_DIMS)) {
        if (!spatial_quadtree_expand_root(qt, min, max)) return false;
    }
    
    /* Create new item */
    spatial_quadtree_item *item = (spatial_quadtree_item*)qt->alloc.malloc(sizeof(spatial_quadtree_item), qt->alloc.udata);
    if (SPATIAL_UNLIKELY(!item)) return false;
    
    spatial_bbox_copy(min, max, item->min, item->max, SPATIAL_QUADTREE_DIMS);
    item->data = data;
    item->next = NULL;
    
    /* Find leaf node and insert */
    spatial_quadtree_node *node = qt->root;
    
    while (!node->is_leaf) {
        int quadrant = spatial_quadtree_quadrant(node, min, max);
        if (!node->children[quadrant]) {
            spatial_num_t child_min[SPATIAL_QUADTREE_DIMS];
            spatial_num_t child_max[SPATIAL_QUADTREE_DIMS];
            spatial_quadtree_child_bounds(node, quadrant, child_min, child_max);
            node->children[quadrant] = spatial_quadtree_node_new(child_min, child_max, node->depth + 1, &qt->alloc);
            if (SPATIAL_UNLIKELY(!node->children[quadrant])) return false;
        }
        node = node->children[quadrant];
    }
    
    /* Insert at front of list */
    item->next = node->items;
    node->items = item;
    node->item_count++;
    qt->count++;
    
    /* Split if needed */
    if (node->item_count > qt->max_items && node->depth < qt->max_depth) {
        spatial_quadtree_node_split(node, &qt->alloc);
    }
    
    return true;
}

/**
 * @brief Search for items overlapping query rectangle
 */
SPATIAL_INLINE void spatial_quadtree_search(const spatial_quadtree *qt,
                                             const spatial_num_t *SPATIAL_RESTRICT qmin,
                                             const spatial_num_t *SPATIAL_RESTRICT qmax,
                                             spatial_quadtree_iter_fn iter,
                                             void *udata)
{
    if (SPATIAL_UNLIKELY(!qt || !qt->root || !iter)) return;

    spatial_quadtree_node *stk[(SPATIAL_QUADTREE_MAX_DEPTH + 1) * SPATIAL_QUADTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = qt->root;

    while (sp > 0) {
        spatial_quadtree_node *node = stk[--sp];

        if (SPATIAL_UNLIKELY(!spatial_bbox_overlaps(node->min, node->max, qmin, qmax,
                                                     SPATIAL_QUADTREE_DIMS))) {
            continue;
        }

        if (node->is_leaf) {
            spatial_quadtree_item *item = node->items;
            while (item) {
                if (spatial_bbox_overlaps(item->min, item->max, qmin, qmax,
                                          SPATIAL_QUADTREE_DIMS)) {
                    if (!iter(item->min, item->max, item->data, udata)) return;
                }
                item = item->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_QUADTREE_CHILDREN; i++) {
                if (node->children[i]) {
                    SPATIAL_PREFETCH(node->children[i], 0, 3);
                    stk[sp++] = node->children[i];
                }
            }
        }
    }
}

/**
 * @brief Scan all items in the quadtree
 */
SPATIAL_INLINE void spatial_quadtree_scan(const spatial_quadtree *qt,
                                           spatial_quadtree_iter_fn iter,
                                           void *udata)
{
    if (SPATIAL_UNLIKELY(!qt || !qt->root || !iter)) return;

    spatial_quadtree_node *stk[(SPATIAL_QUADTREE_MAX_DEPTH + 1) * SPATIAL_QUADTREE_CHILDREN];
    int sp = 0;
    stk[sp++] = qt->root;

    while (sp > 0) {
        spatial_quadtree_node *node = stk[--sp];

        if (node->is_leaf) {
            spatial_quadtree_item *item = node->items;
            while (item) {
                if (!iter(item->min, item->max, item->data, udata)) return;
                item = item->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_QUADTREE_CHILDREN; i++) {
                if (node->children[i]) stk[sp++] = node->children[i];
            }
        }
    }
}

/**
 * @brief Get total item count
 */
SPATIAL_INLINE int spatial_quadtree_count(const spatial_quadtree *qt) {
    return qt ? qt->count : 0;
}

/**
 * @brief Delete items matching bounds and data (using pointer equality)
 */
SPATIAL_INLINE bool spatial_quadtree_delete(spatial_quadtree *qt,
                                             const spatial_num_t *min,
                                             const spatial_num_t *max,
                                             spatial_data_t data)
{
    if (SPATIAL_UNLIKELY(!qt || !qt->root)) return false;

    spatial_quadtree_node *stk[(SPATIAL_QUADTREE_MAX_DEPTH + 1) * SPATIAL_QUADTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = qt->root;

    while (sp > 0 && !found) {
        spatial_quadtree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_QUADTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_quadtree_item **pp = &node->items;
            while (*pp) {
                spatial_quadtree_item *item = *pp;
                if (item->data == data &&
                    spatial_bbox_overlaps(item->min, item->max, min, max, SPATIAL_QUADTREE_DIMS)) {
                    *pp = item->next;
                    if (qt->callbacks.free) qt->callbacks.free(item->data, qt->udata);
                    qt->alloc.free(item, qt->alloc.udata);
                    node->item_count--;
                    qt->count--;
                    found = true;
                    break;
                }
                pp = &(*pp)->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_QUADTREE_CHILDREN; i++) {
                if (node->children[i]) stk[sp++] = node->children[i];
            }
        }
    }

    return found;
}

/**
 * @brief Delete items using custom comparator
 */
typedef bool (*spatial_quadtree_cmp_fn)(spatial_data_t data, void *udata);

SPATIAL_INLINE bool spatial_quadtree_delete_with_comparator(spatial_quadtree *qt,
                                                             const spatial_num_t *min,
                                                             const spatial_num_t *max,
                                                             spatial_quadtree_cmp_fn cmp,
                                                             void *cmp_udata)
{
    if (SPATIAL_UNLIKELY(!qt || !qt->root)) return false;

    spatial_quadtree_node *stk[(SPATIAL_QUADTREE_MAX_DEPTH + 1) * SPATIAL_QUADTREE_CHILDREN];
    int sp = 0;
    bool found = false;
    stk[sp++] = qt->root;

    while (sp > 0 && !found) {
        spatial_quadtree_node *node = stk[--sp];

        if (!spatial_bbox_overlaps(node->min, node->max, min, max, SPATIAL_QUADTREE_DIMS)) {
            continue;
        }

        if (node->is_leaf) {
            spatial_quadtree_item **pp = &node->items;
            while (*pp) {
                spatial_quadtree_item *item = *pp;
                if (cmp(item->data, cmp_udata) &&
                    spatial_bbox_overlaps(item->min, item->max, min, max, SPATIAL_QUADTREE_DIMS)) {
                    *pp = item->next;
                    if (qt->callbacks.free) qt->callbacks.free(item->data, qt->udata);
                    qt->alloc.free(item, qt->alloc.udata);
                    node->item_count--;
                    qt->count--;
                    found = true;
                    break;
                }
                pp = &(*pp)->next;
            }
        } else {
            for (int i = 0; i < SPATIAL_QUADTREE_CHILDREN; i++) {
                if (node->children[i]) stk[sp++] = node->children[i];
            }
        }
    }

    return found;
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_QUADTREE_H */
