/**
 * @file   arena_example.c
 * @brief  Demonstrates the built-in spatial_arena_t with libspatial BVH
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * Zero dependencies — no memento required.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define SPATIAL_NUMTYPE double
#include "libspatial/bvh.h"

#define N_ITEMS 10000
#define N_FRAMES 100

static double rand_range(double min, double max) {
    return min + (max - min) * ((double)rand() / (double)RAND_MAX);
}

int main(void) {
    printf("=== Built-in Arena Example ===\n");
    printf("Items per tree: %d\n", N_ITEMS);
    printf("Frames: %d\n\n", N_FRAMES);

    srand(42);

    spatial_arena_t *arena = spatial_arena_new(4 * 1024 * 1024);
    if (!arena) {
        fprintf(stderr, "Failed to create arena\n");
        return 1;
    }

    /* Frame-based simulation: rebuild BVH each frame in the arena */
    clock_t t0 = clock();
    for (int frame = 0; frame < N_FRAMES; frame++) {
        spatial_arena_reset(arena);

        spatial_bvh *bvh = spatial_bvh_new_with_arena(arena);
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
        int count = 0;
        for (int i = 0; i < bvh->prim_count; i++) {
            if (spatial_bbox_overlaps(bvh->prims[i].min, bvh->prims[i].max, qmin, qmax, 3)) {
                count++;
            }
        }
        (void)count;
    }
    clock_t t1 = clock();
    double arena_ms = (double)(t1 - t0) * 1000.0 / CLOCKS_PER_SEC;
    printf("Arena frame simulation: %.3f ms total (%.3f us/frame)\n\n",
           arena_ms, arena_ms * 1000.0 / N_FRAMES);

    /* Save/restore (cheap CoW) demonstration */
    spatial_arena_reset(arena);
    spatial_bvh *bvh = spatial_bvh_new_with_arena(arena);
    for (int i = 0; i < N_ITEMS; i++) {
        double min[3] = { rand_range(-1000.0, 1000.0),
                          rand_range(-1000.0, 1000.0),
                          rand_range(-1000.0, 1000.0) };
        double max[3] = { min[0] + 1.0, min[1] + 1.0, min[2] + 1.0 };
        spatial_bvh_insert(bvh, min, max, (void*)(intptr_t)i);
    }
    spatial_bvh_build(bvh);

    spatial_arena_save_t save = spatial_arena_save(arena);
    printf("Arena state saved.\n");

    /* Mutate the arena (e.g. speculative work) */
    spatial_bvh *bvh2 = spatial_bvh_new_with_arena(arena);
    for (int i = 0; i < 100; i++) {
        double min[3] = {0, 0, 0};
        double max[3] = {1, 1, 1};
        spatial_bvh_insert(bvh2, min, max, (void*)(intptr_t)i);
    }
    spatial_bvh_build(bvh2);
    printf("Arena used after second BVH.\n");

    spatial_arena_restore(arena, &save);
    printf("Arena restored to saved state.\n\n");

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
        spatial_arena_reset(arena);
        spatial_bvh *abvh = spatial_bvh_new_with_arena(arena);
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
    printf("Arena build (10x):  %.3f ms\n", arena_build_ms);
    printf("Arena speedup:      %.2fx\n\n", malloc_ms / arena_build_ms);

    printf("SUCCESS -- built-in arena zero-overhead mode\n");

    spatial_arena_free(arena);
    return 0;
}
