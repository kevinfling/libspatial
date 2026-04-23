/**
 * @file   spatial.h
 * @brief  Core definitions and utilities for libspatial
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * libspatial — Header-only spatial indexes supporting 8 production-grade structures.
 */

#ifndef SPATIAL_H
#define SPATIAL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>

/* ============================================================================
 * Version Information
 * ============================================================================ */

#define SPATIAL_VERSION_MAJOR 0
#define SPATIAL_VERSION_MINOR 1
#define SPATIAL_VERSION_PATCH 0
#define SPATIAL_VERSION_STRING "0.1.0"

/* ============================================================================
 * Compiler Attributes and Intrinsics
 * ============================================================================ */

#if defined(__GNUC__) || defined(__clang__)
    #define SPATIAL_LIKELY(x)   __builtin_expect(!!(x), 1)
    #define SPATIAL_UNLIKELY(x) __builtin_expect(!!(x), 0)
    #define SPATIAL_INLINE      static inline
    #define SPATIAL_NOINLINE    __attribute__((noinline))
    #define SPATIAL_PURE        __attribute__((pure))
    #define SPATIAL_CONST       __attribute__((const))
    #define SPATIAL_COLD        __attribute__((cold))
    #define SPATIAL_HOT         __attribute__((hot))
    #define SPATIAL_UNUSED      __attribute__((unused))
    #define SPATIAL_ALIGNED(x)  __attribute__((aligned(x)))
    #define SPATIAL_RESTRICT    __restrict__
    #define SPATIAL_PREFETCH(p,rw,loc) __builtin_prefetch((p),(rw),(loc))
#elif defined(_MSC_VER)
    #define SPATIAL_LIKELY(x)   (x)
    #define SPATIAL_UNLIKELY(x) (x)
    #define SPATIAL_INLINE      static inline
    #define SPATIAL_NOINLINE    __declspec(noinline)
    #define SPATIAL_PURE
    #define SPATIAL_CONST
    #define SPATIAL_COLD
    #define SPATIAL_HOT
    #define SPATIAL_UNUSED      __pragma(warning(suppress:4100))
    #define SPATIAL_ALIGNED(x)  __declspec(align(x))
    #define SPATIAL_RESTRICT    __restrict
    #define SPATIAL_PREFETCH(p,rw,loc) _mm_prefetch((const char*)(p),(loc))
#else
    #define SPATIAL_LIKELY(x)   (x)
    #define SPATIAL_UNLIKELY(x) (x)
    #define SPATIAL_INLINE      static inline
    #define SPATIAL_NOINLINE
    #define SPATIAL_PURE
    #define SPATIAL_CONST
    #define SPATIAL_COLD
    #define SPATIAL_HOT
    #define SPATIAL_UNUSED
    #define SPATIAL_ALIGNED(x)
    #define SPATIAL_RESTRICT
    #define SPATIAL_PREFETCH(p,rw,loc) ((void)0)
#endif

/* ============================================================================
 * Configuration Macros with Static Asserts
 * ============================================================================ */

#ifndef SPATIAL_NUMTYPE
    #define SPATIAL_NUMTYPE double
#endif

#ifndef SPATIAL_DIMS
    #define SPATIAL_DIMS 3
#endif

#ifndef SPATIAL_MAXITEMS
    #define SPATIAL_MAXITEMS 64
#endif

#ifndef SPATIAL_DATATYPE
    #define SPATIAL_DATATYPE void*
#endif

/* Static assert for configuration sanity */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
    #define SPATIAL_STATIC_ASSERT(name, cond) _Static_assert((cond), #name)
#else
    #define SPATIAL_STATIC_ASSERT(name, cond) \
        typedef char spatial_static_assert_##name[(cond) ? 1 : -1] SPATIAL_UNUSED
#endif

SPATIAL_STATIC_ASSERT(dims_positive, SPATIAL_DIMS > 0 && SPATIAL_DIMS <= 16);
SPATIAL_STATIC_ASSERT(maxitems_positive, SPATIAL_MAXITEMS > 0 && SPATIAL_MAXITEMS <= 4096);

typedef SPATIAL_NUMTYPE spatial_num_t;
typedef SPATIAL_DATATYPE spatial_data_t;

#include <float.h>

#ifdef __FINITE_MATH_ONLY__
    #define SPATIAL_INFINITY ((spatial_num_t)(sizeof(spatial_num_t) == sizeof(float) ? FLT_MAX : DBL_MAX))
    #define SPATIAL_NEG_INFINITY ((spatial_num_t)(sizeof(spatial_num_t) == sizeof(float) ? -FLT_MAX : -DBL_MAX))
#else
    #define SPATIAL_INFINITY ((spatial_num_t)INFINITY)
    #define SPATIAL_NEG_INFINITY ((spatial_num_t)(-INFINITY))
#endif
#define SPATIAL_EPSILON ((spatial_num_t)1e-9)

/* Alignment for cache lines and SIMD */
#define SPATIAL_CACHE_LINE 64
#define SPATIAL_SIMD_ALIGN 32

/* ============================================================================
 * Atomic Operations
 * ============================================================================ */

#if defined(__STDC_NO_ATOMICS__) || __STDC_VERSION__ < 201112L
    typedef long spatial_atomic_long;
    #define spatial_atomic_load(ptr, order)      (*(ptr))
    #define spatial_atomic_store(ptr, val, order) (*(ptr) = (val))
    #define spatial_atomic_fetch_add(ptr, val, order) ((*(ptr) += (val)) - (val))
    #define spatial_atomic_fetch_sub(ptr, val, order) ((*(ptr) -= (val)) + (val))
    #define SPATIAL_ATOMIC_RELAXED 0
    #define SPATIAL_ATOMIC_SEQ_CST 5
#else
    #include <stdatomic.h>
    typedef _Atomic long spatial_atomic_long;
    #define spatial_atomic_load(ptr, order)      atomic_load_explicit(ptr, order)
    #define spatial_atomic_store(ptr, val, order) atomic_store_explicit(ptr, val, order)
    #define spatial_atomic_fetch_add(ptr, val, order) atomic_fetch_add_explicit(ptr, val, order)
    #define spatial_atomic_fetch_sub(ptr, val, order) atomic_fetch_sub_explicit(ptr, val, order)
    #define SPATIAL_ATOMIC_RELAXED memory_order_relaxed
    #define SPATIAL_ATOMIC_SEQ_CST memory_order_seq_cst
#endif

/* ============================================================================
 * Memory Allocation Interface (with udata for custom allocators like memento)
 * ============================================================================ */

typedef void* (*spatial_malloc_fn)(size_t size, void *udata);
typedef void  (*spatial_free_fn)(void *ptr, void *udata);

typedef struct spatial_allocator {
    spatial_malloc_fn malloc;
    spatial_free_fn   free;
    void *udata;
} spatial_allocator;

SPATIAL_INLINE void* spatial_default_malloc(size_t size, void *udata) {
    (void)udata;
    return malloc(size);
}

SPATIAL_INLINE void spatial_default_free(void *ptr, void *udata) {
    (void)udata;
    free(ptr);
}

SPATIAL_INLINE const spatial_allocator* spatial_default_allocator(void) {
    static const spatial_allocator def = { spatial_default_malloc, spatial_default_free, NULL };
    return &def;
}

/* ============================================================================
 * Built-in Arena (bump allocator)
 * ============================================================================ */

#ifndef SPATIAL_ARENA_BLOCK_SIZE
    #define SPATIAL_ARENA_BLOCK_SIZE ((size_t)65536)
#endif

#ifndef SPATIAL_ARENA_MALLOC
    #define SPATIAL_ARENA_MALLOC(size) malloc(size)
#endif

#ifndef SPATIAL_ARENA_FREE
    #define SPATIAL_ARENA_FREE(ptr) free(ptr)
#endif

typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_arena_block {
    struct spatial_arena_block *next;
    size_t size;
    size_t used;
} spatial_arena_block;

typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_arena {
    spatial_arena_block *head;
} spatial_arena_t;

typedef struct spatial_arena_save {
    spatial_arena_block *block;
    size_t used;
} spatial_arena_save_t;

SPATIAL_INLINE spatial_arena_t* spatial_arena_new(size_t block_size) {
    if (block_size == 0) block_size = SPATIAL_ARENA_BLOCK_SIZE;
    spatial_arena_t *arena = (spatial_arena_t*)SPATIAL_ARENA_MALLOC(sizeof(spatial_arena_t));
    if (SPATIAL_UNLIKELY(!arena)) return NULL;
    arena->head = (spatial_arena_block*)SPATIAL_ARENA_MALLOC(sizeof(spatial_arena_block) + block_size);
    if (SPATIAL_UNLIKELY(!arena->head)) {
        SPATIAL_ARENA_FREE(arena);
        return NULL;
    }
    arena->head->next = NULL;
    arena->head->size = block_size;
    arena->head->used = 0;
    return arena;
}

SPATIAL_INLINE void spatial_arena_free(spatial_arena_t *arena) {
    if (SPATIAL_UNLIKELY(!arena)) return;
    spatial_arena_block *block = arena->head;
    while (block) {
        spatial_arena_block *next = block->next;
        SPATIAL_ARENA_FREE(block);
        block = next;
    }
    SPATIAL_ARENA_FREE(arena);
}

SPATIAL_INLINE void spatial_arena_reset(spatial_arena_t *arena) {
    if (SPATIAL_UNLIKELY(!arena)) return;
    for (spatial_arena_block *block = arena->head; block; block = block->next) {
        block->used = 0;
    }
}

SPATIAL_INLINE void* spatial_arena_alloc(spatial_arena_t *arena, size_t size) {
    if (SPATIAL_UNLIKELY(!arena || size == 0)) return NULL;
    size_t align = SPATIAL_CACHE_LINE;
    size_t pad = (align - (size % align)) % align;
    size_t total = size + pad;

    spatial_arena_block *block = arena->head;
    spatial_arena_block *prev = NULL;
    while (block) {
        if (block->used + total <= block->size) {
            void *ptr = (char*)(block + 1) + block->used;
            block->used += total;
            return ptr;
        }
        prev = block;
        block = block->next;
    }

    size_t new_block_size = SPATIAL_ARENA_BLOCK_SIZE;
    if (total > new_block_size) new_block_size = total;
    spatial_arena_block *new_block = (spatial_arena_block*)SPATIAL_ARENA_MALLOC(sizeof(spatial_arena_block) + new_block_size);
    if (SPATIAL_UNLIKELY(!new_block)) return NULL;
    new_block->next = NULL;
    new_block->size = new_block_size;
    new_block->used = total;
    if (prev) prev->next = new_block;
    else arena->head = new_block;
    return (char*)(new_block + 1);
}

SPATIAL_INLINE spatial_arena_save_t spatial_arena_save(spatial_arena_t *arena) {
    if (!arena || !arena->head) return (spatial_arena_save_t){ NULL, 0 };
    return (spatial_arena_save_t){ arena->head, arena->head->used };
}

SPATIAL_INLINE void spatial_arena_restore(spatial_arena_t *arena, const spatial_arena_save_t *save) {
    if (!arena || !save) {
        if (arena) spatial_arena_reset(arena);
        return;
    }
    if (!save->block) {
        spatial_arena_reset(arena);
        return;
    }
    spatial_arena_block *block = arena->head;
    while (block) {
        if (block == save->block) {
            block->used = save->used;
            spatial_arena_block *to_free = block->next;
            block->next = NULL;
            while (to_free) {
                spatial_arena_block *next = to_free->next;
                SPATIAL_ARENA_FREE(to_free);
                to_free = next;
            }
            return;
        }
        block = block->next;
    }
    /* Save block not found; fall back to full reset */
    spatial_arena_reset(arena);
}

SPATIAL_INLINE void* spatial_arena_malloc(size_t size, void *udata) {
    return spatial_arena_alloc((spatial_arena_t*)udata, size);
}

SPATIAL_INLINE void spatial_arena_adapter_free(void *ptr, void *udata) {
    (void)ptr; (void)udata; /* Arena allocations are bulk-freed */
}

SPATIAL_INLINE spatial_allocator spatial_arena_allocator(spatial_arena_t *arena) {
    spatial_allocator alloc = { spatial_arena_malloc, spatial_arena_adapter_free, arena };
    return alloc;
}

/* ============================================================================
 * Built-in Pool (fixed-size object pool)
 * ============================================================================ */

#ifndef SPATIAL_POOL_BLOCK_COUNT
    #define SPATIAL_POOL_BLOCK_COUNT ((size_t)1024)
#endif

#ifndef SPATIAL_POOL_MALLOC
    #define SPATIAL_POOL_MALLOC(size) malloc(size)
#endif

#ifndef SPATIAL_POOL_FREE
    #define SPATIAL_POOL_FREE(ptr) free(ptr)
#endif

typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_pool_node {
    struct spatial_pool_node *next;
} spatial_pool_node;

typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_pool_block {
    struct spatial_pool_block *next;
} spatial_pool_block;

typedef struct SPATIAL_ALIGNED(SPATIAL_CACHE_LINE) spatial_pool {
    spatial_pool_node *free_list;
    spatial_pool_block *blocks;
    size_t obj_size;
} spatial_pool_t;

SPATIAL_INLINE spatial_pool_t* spatial_pool_new(size_t obj_size, size_t per_block) {
    if (obj_size < sizeof(spatial_pool_node)) obj_size = sizeof(spatial_pool_node);
    if (per_block == 0) per_block = SPATIAL_POOL_BLOCK_COUNT;
    spatial_pool_t *pool = (spatial_pool_t*)SPATIAL_POOL_MALLOC(sizeof(spatial_pool_t));
    if (SPATIAL_UNLIKELY(!pool)) return NULL;
    pool->free_list = NULL;
    pool->blocks = NULL;
    pool->obj_size = obj_size;

    size_t block_bytes = sizeof(spatial_pool_block) + obj_size * per_block;
    spatial_pool_block *block = (spatial_pool_block*)SPATIAL_POOL_MALLOC(block_bytes);
    if (SPATIAL_UNLIKELY(!block)) {
        SPATIAL_POOL_FREE(pool);
        return NULL;
    }
    block->next = NULL;
    pool->blocks = block;

    char *data = (char*)(block + 1);
    for (size_t i = 0; i < per_block; i++) {
        spatial_pool_node *node = (spatial_pool_node*)(data + i * obj_size);
        node->next = pool->free_list;
        pool->free_list = node;
    }
    return pool;
}

SPATIAL_INLINE void spatial_pool_free_all(spatial_pool_t *pool) {
    if (!pool) return;
    spatial_pool_block *block = pool->blocks;
    while (block) {
        spatial_pool_block *next = block->next;
        SPATIAL_POOL_FREE(block);
        block = next;
    }
    pool->blocks = NULL;
    pool->free_list = NULL;
}

SPATIAL_INLINE void* spatial_pool_alloc(spatial_pool_t *pool) {
    if (!pool) return NULL;
    if (!pool->free_list) {
        size_t per_block = SPATIAL_POOL_BLOCK_COUNT;
        size_t block_bytes = sizeof(spatial_pool_block) + pool->obj_size * per_block;
        spatial_pool_block *block = (spatial_pool_block*)SPATIAL_POOL_MALLOC(block_bytes);
        if (!block) return NULL;
        block->next = pool->blocks;
        pool->blocks = block;
        char *data = (char*)(block + 1);
        for (size_t i = 0; i < per_block; i++) {
            spatial_pool_node *node = (spatial_pool_node*)(data + i * pool->obj_size);
            node->next = pool->free_list;
            pool->free_list = node;
        }
    }
    spatial_pool_node *node = pool->free_list;
    pool->free_list = node->next;
    return node;
}

SPATIAL_INLINE void spatial_pool_return(spatial_pool_t *pool, void *ptr) {
    if (!pool || !ptr) return;
    spatial_pool_node *node = (spatial_pool_node*)ptr;
    node->next = pool->free_list;
    pool->free_list = node;
}

SPATIAL_INLINE void* spatial_pool_malloc(size_t size, void *udata) {
    (void)size;
    return spatial_pool_alloc((spatial_pool_t*)udata);
}

SPATIAL_INLINE void spatial_pool_free(void *ptr, void *udata) {
    spatial_pool_return((spatial_pool_t*)udata, ptr);
}

SPATIAL_INLINE spatial_allocator spatial_pool_allocator(spatial_pool_t *pool) {
    spatial_allocator alloc = { spatial_pool_malloc, spatial_pool_free, pool };
    return alloc;
}

/* ============================================================================
 * Memento Integration Helpers
 * ============================================================================ */

#ifdef MEMENTO_H

SPATIAL_INLINE void* spatial_memento_arena_malloc(size_t size, void *udata) {
    return memento_arena_alloc((memento_arena_t*)udata, size, _Alignof(max_align_t));
}

SPATIAL_INLINE void spatial_memento_arena_free(void *ptr, void *udata) {
    (void)ptr; (void)udata; /* Arena allocations are bulk-freed */
}

SPATIAL_INLINE void* spatial_memento_pool_malloc(size_t size, void *udata) {
    (void)size;
    return memento_pool_alloc((memento_pool_t*)udata);
}

SPATIAL_INLINE void spatial_memento_pool_free(void *ptr, void *udata) {
    memento_pool_free((memento_pool_t*)udata, ptr);
}

SPATIAL_INLINE spatial_allocator spatial_memento_arena_allocator(memento_arena_t *arena) {
    spatial_allocator alloc = { spatial_memento_arena_malloc, spatial_memento_arena_free, arena };
    return alloc;
}

SPATIAL_INLINE spatial_allocator spatial_memento_pool_allocator(memento_pool_t *pool) {
    spatial_allocator alloc = { spatial_memento_pool_malloc, spatial_memento_pool_free, pool };
    return alloc;
}

#endif /* MEMENTO_H */

/* ============================================================================
 * Item Callbacks for Deep Copy-on-Write
 * ============================================================================ */

typedef void* (*spatial_item_clone_fn)(const void *item, void *udata);
typedef void  (*spatial_item_free_fn)(void *item, void *udata);

typedef struct spatial_item_callbacks {
    spatial_item_clone_fn clone;
    spatial_item_free_fn  free;
} spatial_item_callbacks;

/* ============================================================================
 * Reference Counting Base Structure
 * ============================================================================ */

typedef struct spatial_refcount {
    spatial_atomic_long count;
} spatial_refcount;

SPATIAL_INLINE void spatial_refcount_init(spatial_refcount *ref) {
    spatial_atomic_store(&ref->count, 1, SPATIAL_ATOMIC_RELAXED);
}

SPATIAL_INLINE long spatial_refcount_increment(const spatial_refcount *ref) {
    union { const spatial_atomic_long *in; spatial_atomic_long *out; } pun;
    pun.in = &ref->count;
    return spatial_atomic_fetch_add(pun.out, 1, SPATIAL_ATOMIC_RELAXED) + 1;
}

SPATIAL_INLINE long spatial_refcount_decrement(spatial_refcount *ref) {
    return spatial_atomic_fetch_sub(&ref->count, 1, SPATIAL_ATOMIC_RELAXED) - 1;
}

SPATIAL_INLINE long spatial_refcount_get(const spatial_refcount *ref) {
    return spatial_atomic_load(&ref->count, SPATIAL_ATOMIC_RELAXED);
}

/* ============================================================================
 * Math Utilities
 * ============================================================================ */

SPATIAL_INLINE SPATIAL_PURE spatial_num_t spatial_min(spatial_num_t a, spatial_num_t b) {
    return (a < b) ? a : b;
}

SPATIAL_INLINE SPATIAL_PURE spatial_num_t spatial_max(spatial_num_t a, spatial_num_t b) {
    return (a > b) ? a : b;
}

SPATIAL_INLINE SPATIAL_PURE spatial_num_t spatial_abs(spatial_num_t a) {
    return (a < 0) ? -a : a;
}

SPATIAL_INLINE SPATIAL_PURE bool spatial_fequal(spatial_num_t a, spatial_num_t b) {
    return spatial_abs(a - b) < SPATIAL_EPSILON;
}

SPATIAL_INLINE SPATIAL_PURE bool spatial_isfinite_num(spatial_num_t a) {
#ifdef __FINITE_MATH_ONLY__
    (void)a;
    return true;
#else
    return isfinite(a);
#endif
}

/* ============================================================================
 * Bounding Box Operations (N-D) with NaN/inf guards
 * ============================================================================ */

SPATIAL_INLINE void spatial_bbox_init(spatial_num_t *SPATIAL_RESTRICT min,
                                       spatial_num_t *SPATIAL_RESTRICT max,
                                       int dims) {
    for (int i = 0; i < dims; i++) {
        min[i] = SPATIAL_INFINITY;
        max[i] = SPATIAL_NEG_INFINITY;
    }
}

SPATIAL_INLINE void spatial_bbox_copy(const spatial_num_t *SPATIAL_RESTRICT src_min,
                                       const spatial_num_t *SPATIAL_RESTRICT src_max,
                                       spatial_num_t *SPATIAL_RESTRICT dst_min,
                                       spatial_num_t *SPATIAL_RESTRICT dst_max,
                                       int dims) {
    for (int i = 0; i < dims; i++) {
        dst_min[i] = src_min[i];
        dst_max[i] = src_max[i];
    }
}

SPATIAL_INLINE void spatial_bbox_merge(const spatial_num_t *SPATIAL_RESTRICT a_min,
                                        const spatial_num_t *SPATIAL_RESTRICT a_max,
                                        const spatial_num_t *SPATIAL_RESTRICT b_min,
                                        const spatial_num_t *SPATIAL_RESTRICT b_max,
                                        spatial_num_t *SPATIAL_RESTRICT out_min,
                                        spatial_num_t *SPATIAL_RESTRICT out_max,
                                        int dims) {
    for (int i = 0; i < dims; i++) {
        out_min[i] = spatial_min(a_min[i], b_min[i]);
        out_max[i] = spatial_max(a_max[i], b_max[i]);
    }
}

SPATIAL_INLINE SPATIAL_PURE bool spatial_bbox_overlaps(const spatial_num_t *SPATIAL_RESTRICT a_min,
                                                         const spatial_num_t *SPATIAL_RESTRICT a_max,
                                                         const spatial_num_t *SPATIAL_RESTRICT b_min,
                                                         const spatial_num_t *SPATIAL_RESTRICT b_max,
                                                         int dims) {
    for (int i = 0; i < dims; i++) {
        if (SPATIAL_LIKELY(a_max[i] < b_min[i] || a_min[i] > b_max[i])) {
            return false;
        }
    }
    return true;
}

SPATIAL_INLINE SPATIAL_PURE bool spatial_bbox_contains(const spatial_num_t *SPATIAL_RESTRICT outer_min,
                                                         const spatial_num_t *SPATIAL_RESTRICT outer_max,
                                                         const spatial_num_t *SPATIAL_RESTRICT inner_min,
                                                         const spatial_num_t *SPATIAL_RESTRICT inner_max,
                                                         int dims) {
    for (int i = 0; i < dims; i++) {
        if (inner_min[i] < outer_min[i] || inner_max[i] > outer_max[i]) {
            return false;
        }
    }
    return true;
}

SPATIAL_INLINE SPATIAL_PURE spatial_num_t spatial_bbox_area(const spatial_num_t *SPATIAL_RESTRICT min,
                                                             const spatial_num_t *SPATIAL_RESTRICT max,
                                                             int dims) {
    spatial_num_t area = (spatial_num_t)1.0;
    for (int i = 0; i < dims; i++) {
        spatial_num_t extent = max[i] - min[i];
        if (extent < (spatial_num_t)0.0) extent = (spatial_num_t)0.0;
        area *= extent;
    }
    return area;
}

SPATIAL_INLINE SPATIAL_PURE spatial_num_t spatial_bbox_margin(const spatial_num_t *SPATIAL_RESTRICT min,
                                                               const spatial_num_t *SPATIAL_RESTRICT max,
                                                               int dims) {
    spatial_num_t margin = (spatial_num_t)0.0;
    for (int i = 0; i < dims; i++) {
        margin += (max[i] - min[i]);
    }
    return (spatial_num_t)2.0 * margin;
}

SPATIAL_INLINE SPATIAL_PURE spatial_num_t spatial_bbox_center_dist_sq(const spatial_num_t *SPATIAL_RESTRICT a_min,
                                                                        const spatial_num_t *SPATIAL_RESTRICT a_max,
                                                                        const spatial_num_t *SPATIAL_RESTRICT b_min,
                                                                        const spatial_num_t *SPATIAL_RESTRICT b_max,
                                                                        int dims) {
    spatial_num_t dist_sq = (spatial_num_t)0.0;
    for (int i = 0; i < dims; i++) {
        spatial_num_t ca = (a_min[i] + a_max[i]) * (spatial_num_t)0.5;
        spatial_num_t cb = (b_min[i] + b_max[i]) * (spatial_num_t)0.5;
        spatial_num_t diff = ca - cb;
        dist_sq += diff * diff;
    }
    return dist_sq;
}

SPATIAL_INLINE bool spatial_bbox_valid(const spatial_num_t *SPATIAL_RESTRICT min,
                                        const spatial_num_t *SPATIAL_RESTRICT max,
                                        int dims) {
    for (int i = 0; i < dims; i++) {
        if (!spatial_isfinite_num(min[i]) || !spatial_isfinite_num(max[i])) {
            return false;
        }
        if (min[i] > max[i]) {
            return false;
        }
    }
    return true;
}

/* ============================================================================
 * Hilbert Curve Utilities
 * ============================================================================ */

SPATIAL_INLINE uint32_t spatial_hilbert_rot(int n, uint32_t x, uint32_t y, uint32_t rx, uint32_t ry) {
    if (ry == 0) {
        if (rx == 1) {
            x = (uint32_t)((1 << n) - 1 - x);
            y = (uint32_t)((1 << n) - 1 - y);
        }
        uint32_t t = x;
        x = y;
        y = t;
    }
    return (x << n) | y;
}

SPATIAL_INLINE uint32_t spatial_hilbert_xy2d(int n, uint32_t x, uint32_t y) {
    uint32_t d = 0;
    for (int s = n - 1; s >= 0; s--) {
        uint32_t rx = (x >> s) & 1;
        uint32_t ry = (y >> s) & 1;
        d += ((uint32_t)1 << (2 * s)) * ((3 * rx) ^ ry);
        x &= ((uint32_t)1 << s) - 1;
        y &= ((uint32_t)1 << s) - 1;
        if (ry == 0) {
            if (rx == 1) {
                x = ((uint32_t)1 << s) - 1 - x;
                y = ((uint32_t)1 << s) - 1 - y;
            }
            uint32_t t = x; x = y; y = t;
        }
    }
    return d;
}

SPATIAL_INLINE uint32_t spatial_hilbert_value_2d(spatial_num_t x, spatial_num_t y,
                                                   spatial_num_t min_x, spatial_num_t min_y,
                                                   spatial_num_t max_x, spatial_num_t max_y) {
    spatial_num_t range_x = max_x - min_x;
    spatial_num_t range_y = max_y - min_y;
    if (range_x < SPATIAL_EPSILON) range_x = SPATIAL_EPSILON;
    if (range_y < SPATIAL_EPSILON) range_y = SPATIAL_EPSILON;
    
    uint32_t ix = (uint32_t)(((x - min_x) / range_x) * 65535.0);
    uint32_t iy = (uint32_t)(((y - min_y) / range_y) * 65535.0);
    return spatial_hilbert_xy2d(16, ix, iy);
}

/* ============================================================================
 * SAH (Surface Area Heuristic) Utilities
 * ============================================================================ */

SPATIAL_INLINE SPATIAL_PURE spatial_num_t spatial_sah_cost(spatial_num_t parent_area,
                                                            spatial_num_t left_area,
                                                            spatial_num_t right_area,
                                                            int left_count,
                                                            int right_count,
                                                            spatial_num_t traversal_cost,
                                                            spatial_num_t intersection_cost) {
    if (parent_area < SPATIAL_EPSILON) return SPATIAL_INFINITY;
    spatial_num_t left_prob = left_area / parent_area;
    spatial_num_t right_prob = right_area / parent_area;
    return traversal_cost + intersection_cost * 
           (left_prob * (spatial_num_t)left_count + right_prob * (spatial_num_t)right_count);
}

/* ============================================================================
 * Quickselect / Median Utilities
 * ============================================================================ */

SPATIAL_INLINE int spatial_partition_double(double *arr, int *indices, int left, int right, int pivot_idx) {
    double pivot_val = arr[pivot_idx];
    int tmp_idx = indices[pivot_idx];
    double tmp_val = arr[pivot_idx];
    
    arr[pivot_idx] = arr[right];
    indices[pivot_idx] = indices[right];
    arr[right] = tmp_val;
    indices[right] = tmp_idx;
    
    int store_idx = left;
    for (int i = left; i < right; i++) {
        if (arr[i] < pivot_val) {
            tmp_val = arr[store_idx];
            tmp_idx = indices[store_idx];
            arr[store_idx] = arr[i];
            indices[store_idx] = indices[i];
            arr[i] = tmp_val;
            indices[i] = tmp_idx;
            store_idx++;
        }
    }
    tmp_val = arr[right];
    tmp_idx = indices[right];
    arr[right] = arr[store_idx];
    indices[right] = indices[store_idx];
    arr[store_idx] = tmp_val;
    indices[store_idx] = tmp_idx;
    return store_idx;
}

SPATIAL_INLINE double spatial_quickselect_double(double *arr, int *indices, int left, int right, int k) {
    while (left < right) {
        int pivot_idx = left + (right - left) / 2;
        pivot_idx = spatial_partition_double(arr, indices, left, right, pivot_idx);
        if (pivot_idx == k) {
            return arr[k];
        } else if (k < pivot_idx) {
            right = pivot_idx - 1;
        } else {
            left = pivot_idx + 1;
        }
    }
    return arr[left];
}

/* ============================================================================
 * Stack Utilities (for non-recursive traversal)
 * ============================================================================ */

typedef struct spatial_stack {
    void **items;
    int    capacity;
    int    top;
    spatial_allocator alloc;
} spatial_stack;

SPATIAL_INLINE spatial_stack* spatial_stack_new(int capacity, const spatial_allocator *alloc) {
    if (!alloc) alloc = spatial_default_allocator();
    spatial_stack *s = (spatial_stack*)alloc->malloc(sizeof(spatial_stack), alloc->udata);
    if (!s) return NULL;
    s->items = (void**)alloc->malloc(sizeof(void*) * (size_t)capacity, alloc->udata);
    if (!s->items) {
        alloc->free(s, alloc->udata);
        return NULL;
    }
    s->capacity = capacity;
    s->top = -1;
    s->alloc = *alloc;
    return s;
}

SPATIAL_INLINE void spatial_stack_free(spatial_stack *s) {
    if (!s) return;
    s->alloc.free(s->items, s->alloc.udata);
    s->alloc.free(s, s->alloc.udata);
}

SPATIAL_INLINE bool spatial_stack_push(spatial_stack *s, void *item) {
    if (SPATIAL_UNLIKELY(s->top + 1 >= s->capacity)) return false;
    s->items[++s->top] = item;
    return true;
}

SPATIAL_INLINE void* spatial_stack_pop(spatial_stack *s) {
    if (SPATIAL_UNLIKELY(s->top < 0)) return NULL;
    return s->items[s->top--];
}

SPATIAL_INLINE bool spatial_stack_empty(const spatial_stack *s) {
    return s->top < 0;
}

/* ============================================================================
 * Priority Queue (for k-NN searches)
 * ============================================================================ */

typedef struct spatial_pq_node {
    void          *item;
    spatial_num_t  dist;
} spatial_pq_node;

typedef struct spatial_priority_queue {
    spatial_pq_node  *nodes;
    int               size;
    int               capacity;
    spatial_allocator alloc;
} spatial_priority_queue;

SPATIAL_INLINE spatial_priority_queue* spatial_pq_new(int capacity, const spatial_allocator *alloc) {
    if (!alloc) alloc = spatial_default_allocator();
    spatial_priority_queue *pq = (spatial_priority_queue*)alloc->malloc(sizeof(spatial_priority_queue), alloc->udata);
    if (!pq) return NULL;
    pq->nodes = (spatial_pq_node*)alloc->malloc(sizeof(spatial_pq_node) * (size_t)(capacity + 1), alloc->udata);
    if (!pq->nodes) {
        alloc->free(pq, alloc->udata);
        return NULL;
    }
    pq->size = 0;
    pq->capacity = capacity;
    pq->alloc = *alloc;
    return pq;
}

SPATIAL_INLINE void spatial_pq_free(spatial_priority_queue *pq) {
    if (!pq) return;
    pq->alloc.free(pq->nodes, pq->alloc.udata);
    pq->alloc.free(pq, pq->alloc.udata);
}

SPATIAL_INLINE void spatial_pq_swap(spatial_pq_node *a, spatial_pq_node *b) {
    spatial_pq_node t = *a; *a = *b; *b = t;
}

SPATIAL_INLINE void spatial_pq_heapify_up(spatial_priority_queue *pq, int idx) {
    while (idx > 1 && pq->nodes[idx].dist < pq->nodes[idx / 2].dist) {
        spatial_pq_swap(&pq->nodes[idx], &pq->nodes[idx / 2]);
        idx /= 2;
    }
}

SPATIAL_INLINE void spatial_pq_heapify_down(spatial_priority_queue *pq, int idx) {
    int size = pq->size;
    while (2 * idx <= size) {
        int j = 2 * idx;
        if (j < size && pq->nodes[j + 1].dist < pq->nodes[j].dist) j++;
        if (pq->nodes[idx].dist <= pq->nodes[j].dist) break;
        spatial_pq_swap(&pq->nodes[idx], &pq->nodes[j]);
        idx = j;
    }
}

SPATIAL_INLINE bool spatial_pq_push(spatial_priority_queue *pq, void *item, spatial_num_t dist) {
    if (pq->size >= pq->capacity) return false;
    pq->size++;
    pq->nodes[pq->size] = (spatial_pq_node){ item, dist };
    spatial_pq_heapify_up(pq, pq->size);
    return true;
}

SPATIAL_INLINE bool spatial_pq_pop(spatial_priority_queue *pq, void **item, spatial_num_t *dist) {
    if (pq->size == 0) return false;
    *item = pq->nodes[1].item;
    *dist = pq->nodes[1].dist;
    pq->nodes[1] = pq->nodes[pq->size--];
    spatial_pq_heapify_down(pq, 1);
    return true;
}

SPATIAL_INLINE bool spatial_pq_empty(const spatial_priority_queue *pq) {
    return pq->size == 0;
}

SPATIAL_INLINE spatial_num_t spatial_pq_min_dist(const spatial_priority_queue *pq) {
    return (pq->size > 0) ? pq->nodes[1].dist : SPATIAL_INFINITY;
}

#ifdef __cplusplus
}
#endif

#endif /* SPATIAL_H */
