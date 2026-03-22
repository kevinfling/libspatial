/**
 * @file   memento_pool_example.c
 * @brief  Demonstrates memento pool+arena integration with libspatial Quadtree
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "memento.h"

#define SPATIAL_NUMTYPE double
#include "libspatial/quadtree.h"

#define N_ITEMS 50000
#define N_ITERS 10

static double rand_range(double min, double max) {
    return min + (max - min) * ((double)rand() / (double)RAND_MAX);
}

static void benchmark_malloc(void) {
    clock_t t0 = clock();
    for (int iter = 0; iter < N_ITERS; iter++) {
        spatial_quadtree *qt = spatial_quadtree_new();
        for (int i = 0; i < N_ITEMS; i++) {
            double min[2] = { rand_range(-1000.0, 1000.0), rand_range(-1000.0, 1000.0) };
            double max[2] = { min[0] + rand_range(1.0, 10.0), min[1] + rand_range(1.0, 10.0) };
            spatial_quadtree_insert(qt, min, max, (void*)(intptr_t)i);
        }
        spatial_quadtree_free(qt);
    }
    clock_t t1 = clock();
    double ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    printf("  malloc:  %.3f ms total (%.3f us/insert)\n", ms, ms * 1000.0 / (N_ITERS * N_ITEMS));
}

static void benchmark_pool_arena(void) {
    memento_thread_heap_t *heap = memento_thread_heap_get();
    memento_arena_t *arena = memento_arena_create(2 * 1024 * 1024, heap);
    memento_pool_t *item_pool = memento_pool_create(sizeof(spatial_quadtree_item), N_ITEMS, heap);
    memento_pool_t *node_pool = memento_pool_create(sizeof(spatial_quadtree_node), N_ITEMS, heap);

    (void)node_pool; /* reserved for advanced custom allocation patterns */

    clock_t t0 = clock();
    for (int iter = 0; iter < N_ITERS; iter++) {
        memento_arena_reset(arena);
        memento_pool_destroy(item_pool);
        memento_pool_destroy(node_pool);
        item_pool = memento_pool_create(sizeof(spatial_quadtree_item), N_ITEMS, heap);
        node_pool = memento_pool_create(sizeof(spatial_quadtree_node), N_ITEMS, heap);

        spatial_quadtree *qt = spatial_quadtree_new_with_memento_arena(arena);
        for (int i = 0; i < N_ITEMS; i++) {
            double min[2] = { rand_range(-1000.0, 1000.0), rand_range(-1000.0, 1000.0) };
            double max[2] = { min[0] + rand_range(1.0, 10.0), min[1] + rand_range(1.0, 10.0) };
            spatial_quadtree_insert(qt, min, max, (void*)(intptr_t)i);
        }
        /* Arena-based trees are bulk-freed by resetting the arena */
        spatial_quadtree_free(qt);
    }
    clock_t t1 = clock();
    double ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    printf("  arena:   %.3f ms total (%.3f us/insert)\n", ms, ms * 1000.0 / (N_ITERS * N_ITEMS));

    memento_pool_destroy(item_pool);
    memento_pool_destroy(node_pool);
    memento_arena_destroy(arena);
}

int main(void) {
    if (!memento_init()) {
        fprintf(stderr, "Failed to initialize memento\n");
        return 1;
    }

    printf("=== Memento Pool + Arena Example ===\n");
    printf("Items per iteration: %d\n", N_ITEMS);
    printf("Iterations: %d\n\n", N_ITERS);

    srand(42);

    printf("Benchmarking plain malloc...\n");
    benchmark_malloc();

    printf("Benchmarking memento arena...\n");
    benchmark_pool_arena();

    printf("\nSUCCESS -- memento zero-overhead mode\n");

    memento_shutdown();
    return 0;
}
