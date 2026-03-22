/**
 * @file   memento_arena_example.c
 * @brief  Demonstrates memento arena integration with libspatial BVH
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
#include "libspatial/bvh.h"

#define N_ITEMS 10000
#define N_FRAMES 100

static double rand_range(double min, double max) {
    return min + (max - min) * ((double)rand() / (double)RAND_MAX);
}

int main(void) {
    if (!memento_init()) {
        fprintf(stderr, "Failed to initialize memento\n");
        return 1;
    }

    memento_thread_heap_t *heap = memento_thread_heap_get();
    memento_arena_t *arena = memento_arena_create(4 * 1024 * 1024, heap);
    if (!arena) {
        fprintf(stderr, "Failed to create arena\n");
        return 1;
    }

    printf("=== Memento Arena Example ===\n");
    printf("Initial arena capacity: %zu bytes\n", memento_arena_capacity(arena));
    printf("Items per tree: %d\n", N_ITEMS);
    printf("Frames: %d\n\n", N_FRAMES);

    srand(42);

    /* Frame-based simulation: rebuild BVH each frame in the arena */
    clock_t t0 = clock();
    for (int frame = 0; frame < N_FRAMES; frame++) {
        memento_arena_reset(arena);

        spatial_bvh *bvh = spatial_bvh_new_with_memento_arena(arena);
        if (!bvh) {
            fprintf(stderr, "Failed to create BVH\n");
            return 1;
        }

        for (int i = 0; i < N_ITEMS; i++) {
            double min[3] = { rand_range(-1000.0, 1000.0),
                              rand_range(-1000.0, 1000.0),
                              rand_range(-1000.0, 1000.0) };
            double max[3] = { min[0] + rand_range(1.0, 10.0),
                              min[1] + rand_range(1.0, 10.0),
                              min[2] + rand_range(1.0, 10.0) };
            spatial_bvh_insert(bvh, min, max, (void*)(intptr_t)i);
        }

        if (!spatial_bvh_build(bvh)) {
            fprintf(stderr, "Failed to build BVH\n");
            return 1;
        }

        /* Simulate a search */
        double qmin[3] = {-100.0, -100.0, -100.0};
        double qmax[3] = {100.0, 100.0, 100.0};
        /* Simple search counter */
        struct { int count; } ctx = {0};
        for (int i = 0; i < bvh->prim_count; i++) {
            if (spatial_bbox_overlaps(bvh->prims[i].min, bvh->prims[i].max, qmin, qmax, 3)) {
                ctx.count++;
            }
        }
    }
    clock_t t1 = clock();
    double arena_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    printf("Arena frame simulation: %.3f ms total (%.3f us/frame)\n",
           arena_ms, arena_ms * 1000.0 / N_FRAMES);
    printf("Arena used after frames: %zu bytes\n\n", memento_arena_used(arena));

    /* CoW clone demonstration */
    memento_arena_reset(arena);
    spatial_bvh *bvh = spatial_bvh_new_with_memento_arena(arena);
    for (int i = 0; i < N_ITEMS; i++) {
        double min[3] = { rand_range(-1000.0, 1000.0),
                          rand_range(-1000.0, 1000.0),
                          rand_range(-1000.0, 1000.0) };
        double max[3] = { min[0] + 1.0, min[1] + 1.0, min[2] + 1.0 };
        spatial_bvh_insert(bvh, min, max, (void*)(intptr_t)i);
    }
    spatial_bvh_build(bvh);

    memento_arena_save_t save = memento_arena_save(arena);
    printf("Arena used before clone: %zu bytes\n", memento_arena_used(arena));

    /* "Clone" by saving state, doing more work, then restoring */
    spatial_bvh *bvh2 = spatial_bvh_new_with_memento_arena(arena);
    for (int i = 0; i < 100; i++) {
        double min[3] = {0, 0, 0};
        double max[3] = {1, 1, 1};
        spatial_bvh_insert(bvh2, min, max, (void*)(intptr_t)i);
    }
    spatial_bvh_build(bvh2);
    printf("Arena used after second BVH: %zu bytes\n", memento_arena_used(arena));

    memento_arena_restore(arena, &save);
    printf("Arena used after restore: %zu bytes\n\n", memento_arena_used(arena));

    /* Benchmark: arena vs default malloc */
    t0 = clock();
    for (int iter = 0; iter < 10; iter++) {
        spatial_bvh *mbvh = spatial_bvh_new();
        for (int i = 0; i < N_ITEMS; i++) {
            double min[3] = { rand_range(-1000.0, 1000.0),
                              rand_range(-1000.0, 1000.0),
                              rand_range(-1000.0, 1000.0) };
            double max[3] = { min[0] + 1.0, min[1] + 1.0, min[2] + 1.0 };
            spatial_bvh_insert(mbvh, min, max, (void*)(intptr_t)i);
        }
        spatial_bvh_build(mbvh);
        spatial_bvh_free(mbvh);
    }
    t1 = clock();
    double malloc_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    printf("Malloc build (10x): %.3f ms\n", malloc_ms);

    t0 = clock();
    for (int iter = 0; iter < 10; iter++) {
        memento_arena_reset(arena);
        spatial_bvh *abvh = spatial_bvh_new_with_memento_arena(arena);
        for (int i = 0; i < N_ITEMS; i++) {
            double min[3] = { rand_range(-1000.0, 1000.0),
                              rand_range(-1000.0, 1000.0),
                              rand_range(-1000.0, 1000.0) };
            double max[3] = { min[0] + 1.0, min[1] + 1.0, min[2] + 1.0 };
            spatial_bvh_insert(abvh, min, max, (void*)(intptr_t)i);
        }
        spatial_bvh_build(abvh);
    }
    t1 = clock();
    double arena_build_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    printf("Arena build (10x): %.3f ms\n", arena_build_ms);
    printf("Arena speedup: %.2fx\n\n", malloc_ms / arena_build_ms);

    printf("SUCCESS -- memento zero-overhead mode\n");

    memento_arena_destroy(arena);
    memento_shutdown();
    return 0;
}
