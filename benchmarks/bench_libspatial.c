/**
 * @file   bench_libspatial.c
 * @brief  Performance benchmarks for libspatial
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * Measures:
 * - Insertion performance (items/second)
 * - Search performance (queries/second)
 * - Memory usage (approximate)
 * - Scalability from 1K to 1M items
 * Build with: -O3 -march=native -DNDEBUG
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>

#ifdef MEMENTO_IMPLEMENTATION
#include "memento.h"
#endif

/* Platform-specific timing */
#if defined(__x86_64__) || defined(_M_X64)
    #include <x86intrin.h>
    static inline uint64_t rdtsc(void) {
        uint32_t lo, hi;
        __asm__ volatile("lfence\n\trdtsc" : "=a"(lo), "=d"(hi));
        return ((uint64_t)hi << 32) | lo;
    }
    static inline double timer_hz(void) {
        /* Approximate; rdtsc runs at CPU clock. Calibrate if precision matters. */
        return 3.0e9;
    }
#elif defined(__aarch64__)
    static inline uint64_t rdtsc(void) {
        uint64_t val;
        __asm__ volatile("mrs %0, cntvct_el0" : "=r"(val));
        return val;
    }
    static inline double timer_hz(void) {
        uint64_t freq;
        __asm__ volatile("mrs %0, cntfrq_el0" : "=r"(freq));
        return (double)freq;
    }
#else
    static inline uint64_t rdtsc(void) {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return (uint64_t)ts.tv_sec * 1000000000ULL + ts.tv_nsec;
    }
    static inline double timer_hz(void) {
        return 1.0e9;
    }
#endif

static void rdtsc_warmup(void) {
    for (int i = 0; i < 10000; i++) {
        (void)rdtsc();
    }
}

#include "libspatial/libspatial.h"

/* Benchmark configuration */
#define BENCH_MIN_ITEMS     1000
#define BENCH_MAX_ITEMS     1000000
#define BENCH_NUM_QUERIES   10000
#define BENCH_WARMUP        100

typedef struct {
    const char *name;
    uint64_t insert_cycles;
    uint64_t search_cycles;
    uint64_t scan_cycles;
    double insert_ns;
    double search_ns;
    double scan_ns;
} bench_result;

static double cycles_to_ns(uint64_t cycles) {
    return (double)cycles / (timer_hz() / 1.0e9);
}

static void print_header(void) {
    printf("\n");
    printf("%-20s %12s %12s %12s %12s %12s %12s\n",
           "Structure", "Insert(cyc)", "Search(cyc)", "Scan(cyc)",
           "Insert(ns)", "Search(ns)", "Scan(ns)");
    printf("%-20s %12s %12s %12s %12s %12s %12s\n",
           "--------------------", "------------", "------------", "------------",
           "------------", "------------", "------------");
}

static void print_result(const bench_result *r) {
    printf("%-20s %12.1f %12.1f %12.1f %12.1f %12.1f %12.1f\n",
           r->name,
           (double)r->insert_cycles,
           (double)r->search_cycles,
           (double)r->scan_cycles,
           r->insert_ns,
           r->search_ns,
           r->scan_ns);
}

/* Generate random AABB */
static void random_aabb(double *min, double *max, int dims, double range) {
    for (int i = 0; i < dims; i++) {
        min[i] = ((double)rand() / RAND_MAX) * range;
        double size = ((double)rand() / RAND_MAX) * range * 0.01;
        if (size < 0.1) size = 0.1;
        max[i] = min[i] + size;
    }
}

/* Callback that just counts */
static volatile int g_hit_count = 0;
static bool count_callback(const spatial_num_t *min, const spatial_num_t *max,
                           spatial_data_t data, void *udata) {
    (void)min; (void)max; (void)data; (void)udata;
    g_hit_count++;
    return true;
}

/* ============================================================================
 * Quadtree Benchmark
 * ============================================================================ */

static bench_result bench_quadtree(int num_items) {
    bench_result r = {"Quadtree", 0, 0, 0, 0, 0, 0};
    
    spatial_quadtree *qt = spatial_quadtree_new();
    if (!qt) return r;
    
    double *mins = malloc(num_items * 2 * 2 * sizeof(double));
    double *maxs = mins + num_items * 2;
    
    /* Generate data */
    for (int i = 0; i < num_items; i++) {
        random_aabb(&mins[i * 2], &maxs[i * 2], 2, 10000.0);
    }
    
    /* Warmup */
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_quadtree_insert(qt, &mins[i * 2], &maxs[i * 2], (void*)(intptr_t)i);
    }
    
    spatial_quadtree_free(qt);
    qt = spatial_quadtree_new();
    
    /* Benchmark insert */
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_quadtree_insert(qt, &mins[i * 2], &maxs[i * 2], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    /* Benchmark search */
    double qmin[2], qmax[2];
    random_aabb(qmin, qmax, 2, 10000.0);
    
    g_hit_count = 0;
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_quadtree_search(qt, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    /* Benchmark scan */
    g_hit_count = 0;
    start = rdtsc();
    spatial_quadtree_scan(qt, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(mins);
    spatial_quadtree_free(qt);
    
    return r;
}

/* ============================================================================
 * Octree Benchmark
 * ============================================================================ */

static bench_result bench_octree(int num_items) {
    bench_result r = {"Octree", 0, 0, 0, 0, 0, 0};
    
    spatial_octree *ot = spatial_octree_new();
    if (!ot) return r;
    
    double *mins = malloc(num_items * 3 * 2 * sizeof(double));
    double *maxs = mins + num_items * 3;
    
    for (int i = 0; i < num_items; i++) {
        random_aabb(&mins[i * 3], &maxs[i * 3], 3, 10000.0);
    }
    
    /* Warmup */
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_octree_insert(ot, &mins[i * 3], &maxs[i * 3], (void*)(intptr_t)i);
    }
    
    spatial_octree_free(ot);
    ot = spatial_octree_new();
    
    /* Insert */
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_octree_insert(ot, &mins[i * 3], &maxs[i * 3], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    /* Search */
    double qmin[3], qmax[3];
    random_aabb(qmin, qmax, 3, 10000.0);
    
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_octree_search(ot, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    /* Scan */
    start = rdtsc();
    spatial_octree_scan(ot, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(mins);
    spatial_octree_free(ot);
    
    return r;
}

/* ============================================================================
 * Orthtree Benchmark (4D default)
 * ============================================================================ */

static bench_result bench_orthtree(int num_items) {
    bench_result r = {"Orthtree", 0, 0, 0, 0, 0, 0};
    
    spatial_orthtree *ot = spatial_orthtree_new();
    if (!ot) return r;
    
    double *mins = malloc(num_items * 4 * 2 * sizeof(double));
    double *maxs = mins + num_items * 4;
    
    for (int i = 0; i < num_items; i++) {
        random_aabb(&mins[i * 4], &maxs[i * 4], 4, 10000.0);
    }
    
    /* Warmup */
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_orthtree_insert(ot, &mins[i * 4], &maxs[i * 4], (void*)(intptr_t)i);
    }
    spatial_orthtree_free(ot);
    ot = spatial_orthtree_new();
    
    /* Insert */
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_orthtree_insert(ot, &mins[i * 4], &maxs[i * 4], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    /* Search */
    double qmin[4], qmax[4];
    random_aabb(qmin, qmax, 4, 10000.0);
    
    g_hit_count = 0;
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_orthtree_search(ot, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    /* Scan */
    start = rdtsc();
    spatial_orthtree_scan(ot, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(mins);
    spatial_orthtree_free(ot);
    
    return r;
}

/* ============================================================================
 * KD-Tree Benchmark
 * ============================================================================ */

static bench_result bench_kdtree(int num_items) {
    bench_result r = {"KD-Tree", 0, 0, 0, 0, 0, 0};
    
    spatial_kdtree *kt = spatial_kdtree_new();
    if (!kt) return r;
    
    double *points = malloc(num_items * 3 * sizeof(double));
    
    for (int i = 0; i < num_items; i++) {
        random_aabb(&points[i * 3], &points[i * 3], 3, 10000.0);
    }
    
    /* Warmup */
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_kdtree_insert(kt, &points[i * 3], &points[i * 3], (void*)(intptr_t)i);
    }
    spatial_kdtree_free(kt);
    kt = spatial_kdtree_new();
    
    /* Insert */
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_kdtree_insert(kt, &points[i * 3], &points[i * 3], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    /* Search */
    double qmin[3], qmax[3];
    random_aabb(qmin, qmax, 3, 10000.0);
    
    g_hit_count = 0;
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_kdtree_search(kt, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    /* Scan */
    start = rdtsc();
    spatial_kdtree_scan(kt, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(points);
    spatial_kdtree_free(kt);
    
    return r;
}

/* ============================================================================
 * VP-Tree Benchmark
 * ============================================================================ */

static bench_result bench_vptree(int num_items) {
    bench_result r = {"VP-Tree", 0, 0, 0, 0, 0, 0};
    
    spatial_vptree *vp = spatial_vptree_new();
    if (!vp) return r;
    
    double *points = malloc(num_items * 3 * sizeof(double));
    
    for (int i = 0; i < num_items; i++) {
        random_aabb(&points[i * 3], &points[i * 3], 3, 10000.0);
    }
    
    /* Warmup */
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_vptree_insert(vp, &points[i * 3], &points[i * 3], (void*)(intptr_t)i);
    }
    spatial_vptree_free(vp);
    vp = spatial_vptree_new();
    
    /* Insert */
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_vptree_insert(vp, &points[i * 3], &points[i * 3], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    /* Search */
    double qmin[3], qmax[3];
    random_aabb(qmin, qmax, 3, 10000.0);
    
    g_hit_count = 0;
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_vptree_search(vp, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    /* Scan */
    start = rdtsc();
    spatial_vptree_scan(vp, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(points);
    spatial_vptree_free(vp);
    
    return r;
}

/* ============================================================================
 * Hilbert R-Tree Benchmark
 * ============================================================================ */

static bench_result bench_hilbertrtree(int num_items) {
    bench_result r = {"Hilbert R-Tree", 0, 0, 0, 0, 0, 0};
    
    spatial_hilbertrtree *rt = spatial_hilbertrtree_new();
    if (!rt) return r;
    
    double *mins = malloc(num_items * 2 * 2 * sizeof(double));
    double *maxs = mins + num_items * 2;
    
    for (int i = 0; i < num_items; i++) {
        random_aabb(&mins[i * 2], &maxs[i * 2], 2, 10000.0);
    }
    
    /* Warmup */
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_hilbertrtree_insert(rt, &mins[i * 2], &maxs[i * 2], (void*)(intptr_t)i);
    }
    spatial_hilbertrtree_free(rt);
    rt = spatial_hilbertrtree_new();
    
    /* Insert */
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_hilbertrtree_insert(rt, &mins[i * 2], &maxs[i * 2], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    spatial_hilbertrtree_rebuild(rt);
    
    /* Search */
    double qmin[2], qmax[2];
    random_aabb(qmin, qmax, 2, 10000.0);
    
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_hilbertrtree_search(rt, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    /* Scan */
    start = rdtsc();
    spatial_hilbertrtree_scan(rt, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(mins);
    spatial_hilbertrtree_free(rt);
    
    return r;
}

/* ============================================================================
 * BSP-Tree Benchmark
 * ============================================================================ */
static bench_result bench_bsptree(int num_items) {
    bench_result r = {"BSP-Tree", 0, 0, 0, 0, 0, 0};
    
    spatial_bsptree *bsp = spatial_bsptree_new();
    if (!bsp) return r;
    
    double *mins = malloc(num_items * 3 * 2 * sizeof(double));
    double *maxs = mins + num_items * 3;
    
    for (int i = 0; i < num_items; i++) {
        random_aabb(&mins[i * 3], &maxs[i * 3], 3, 10000.0);
    }
    
    // Warmup
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_bsptree_insert(bsp, &mins[i * 3], &maxs[i * 3], (void*)(intptr_t)i);
    }
    spatial_bsptree_free(bsp);
    bsp = spatial_bsptree_new();
    
    // Insert 
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_bsptree_insert(bsp, &mins[i * 3], &maxs[i * 3], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    spatial_bsptree_rebuild(bsp);
    
    // Search
    double qmin[3], qmax[3];
    random_aabb(qmin, qmax, 3, 10000.0);
    
    g_hit_count = 0;
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_bsptree_search(bsp, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    // Scan
    start = rdtsc();
    spatial_bsptree_scan(bsp, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(mins);
    spatial_bsptree_free(bsp);

    return r;
}

/* ============================================================================
 * BVH Benchmark
 * ============================================================================ */
static bench_result bench_bvh(int num_items) {
    bench_result r = {"BVH", 0, 0, 0, 0, 0, 0};
    
    spatial_bvh *bvh = spatial_bvh_new();
    if (!bvh) return r;
    
    double *mins = malloc(num_items * 3 * 2 * sizeof(double));
    double *maxs = mins + num_items * 3;
    
    for (int i = 0; i < num_items; i++) {
        random_aabb(&mins[i * 3], &maxs[i * 3], 3, 10000.0);
    }
    
    // Warmup
    for (int i = 0; i < BENCH_WARMUP && i < num_items; i++) {
        spatial_bvh_insert(bvh, &mins[i * 3], &maxs[i * 3], (void*)(intptr_t)i);
    }
    spatial_bvh_free(bvh);
    bvh = spatial_bvh_new();
    
    // Insert
    uint64_t start = rdtsc();
    for (int i = 0; i < num_items; i++) {
        spatial_bvh_insert(bvh, &mins[i * 3], &maxs[i * 3], (void*)(intptr_t)i);
    }
    uint64_t end = rdtsc();
    r.insert_cycles = (end - start) / num_items;
    r.insert_ns = cycles_to_ns(r.insert_cycles);
    
    // Build
    spatial_bvh_build(bvh);
    
    // Search
    double qmin[3], qmax[3];
    random_aabb(qmin, qmax, 3, 10000.0);
    
    g_hit_count = 0;
    start = rdtsc();
    for (int i = 0; i < BENCH_NUM_QUERIES; i++) {
        spatial_bvh_search(bvh, qmin, qmax, count_callback, NULL);
    }
    end = rdtsc();
    r.search_cycles = (end - start) / BENCH_NUM_QUERIES;
    r.search_ns = cycles_to_ns(r.search_cycles);
    
    // Scan
    start = rdtsc();
    spatial_bvh_scan(bvh, count_callback, NULL);
    end = rdtsc();
    r.scan_cycles = (end - start) / num_items;
    r.scan_ns = cycles_to_ns(r.scan_cycles);
    
    free(mins);
    spatial_bvh_free(bvh);

    return r;
}
/* ============================================================================
 * Main
 * ============================================================================ */

int main(int argc, char **argv) {
    (void)argc; (void)argv;
    
    srand(42);
    rdtsc_warmup();
    
    int num_items = 1000000;
    if (argc > 1) {
        num_items = atoi(argv[1]);
        if (num_items < BENCH_MIN_ITEMS) num_items = BENCH_MIN_ITEMS;
        if (num_items > BENCH_MAX_ITEMS) num_items = BENCH_MAX_ITEMS;
    }
    
    printf("\n");
    printf("=================================================================\n");
    printf("libspatial Performance Benchmark v%s\n", LIBSPATIAL_VERSION_STRING);
    printf("=================================================================\n");
    printf("Items: %d, Queries: %d\n", num_items, BENCH_NUM_QUERIES);
    printf("(Timer ticks; ns conversion uses actual counter frequency)\n");
    
    print_header();
    
    bench_result r;
    
    r = bench_quadtree(num_items);
    print_result(&r);
    
    r = bench_octree(num_items);
    print_result(&r);
    
    r = bench_orthtree(num_items);
    print_result(&r);
    
    r = bench_kdtree(num_items);
    print_result(&r);
    
    r = bench_vptree(num_items);
    print_result(&r);
    
    r = bench_hilbertrtree(num_items);
    print_result(&r);
    r = bench_bsptree(num_items);
    print_result(&r);

    r = bench_bvh(num_items);
    print_result(&r);
    printf("\n");
    printf("=================================================================\n");
    printf("Target: < 200 ns/op search on 1M items (AVX-ready layout)\n");
    printf("=================================================================\n\n");

    /* Allocator comparison */
    {
        const int alloc_items = num_items < 100000 ? num_items : 100000;
        double *amins = (double*)malloc(sizeof(double) * alloc_items * 3);
        double *amaxs = (double*)malloc(sizeof(double) * alloc_items * 3);
        if (amins && amaxs) {
            for (int i = 0; i < alloc_items; i++) {
                random_aabb(&amins[i * 3], &amaxs[i * 3], 3, 10000.0);
            }
            const int iters = 10;

            uint64_t start, end;

            /* malloc */
            start = rdtsc();
            for (int iter = 0; iter < iters; iter++) {
                spatial_bvh *bvh = spatial_bvh_new();
                for (int i = 0; i < alloc_items; i++) {
                    spatial_bvh_insert(bvh, &amins[i * 3], &amaxs[i * 3], (void*)(intptr_t)i);
                }
                spatial_bvh_build(bvh);
                spatial_bvh_free(bvh);
            }
            end = rdtsc();
            double malloc_ms = cycles_to_ns(end - start) / 1e6;

            /* built-in arena */
            spatial_arena_t *arena = spatial_arena_new(4 * 1024 * 1024);
            start = rdtsc();
            for (int iter = 0; iter < iters; iter++) {
                spatial_arena_reset(arena);
                spatial_bvh *bvh = spatial_bvh_new_with_arena(arena);
                for (int i = 0; i < alloc_items; i++) {
                    spatial_bvh_insert(bvh, &amins[i * 3], &amaxs[i * 3], (void*)(intptr_t)i);
                }
                spatial_bvh_build(bvh);
            }
            end = rdtsc();
            double arena_ms = cycles_to_ns(end - start) / 1e6;

#ifdef MEMENTO_H
            /* memento arena */
            memento_init();
            memento_thread_heap_t *heap = memento_thread_heap_get();
            memento_arena_t *marena = memento_arena_create(4 * 1024 * 1024, heap);
            start = rdtsc();
            for (int iter = 0; iter < iters; iter++) {
                memento_arena_reset(marena);
                spatial_bvh *bvh = spatial_bvh_new_with_memento_arena(marena);
                for (int i = 0; i < alloc_items; i++) {
                    spatial_bvh_insert(bvh, &amins[i * 3], &amaxs[i * 3], (void*)(intptr_t)i);
                }
                spatial_bvh_build(bvh);
            }
            end = rdtsc();
            double memento_ms = cycles_to_ns(end - start) / 1e6;
            memento_arena_destroy(marena);
            memento_shutdown();
#endif

            printf("=== Allocator Comparison (%d items, %d iters) ===\n", alloc_items, iters);
            printf("  malloc:        %.3f ms\n", malloc_ms);
            printf("  spatial_arena: %.3f ms (%.2fx vs malloc)\n", arena_ms, malloc_ms / arena_ms);
#ifdef MEMENTO_H
            printf("  memento_arena: %.3f ms (%.2fx vs malloc)\n", memento_ms, malloc_ms / memento_ms);
#endif
            printf("\n");

            spatial_arena_free(arena);
        }
        free(amins);
        free(amaxs);
    }

    return 0;
}
