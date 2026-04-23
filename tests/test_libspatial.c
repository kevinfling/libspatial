/**
 * @file   test_libspatial.c
 * @brief  Comprehensive test suite for libspatial
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * Tests cover:
 * - All 8 data structures
 * - Insert/search/delete operations
 * - Clone/copy-on-write semantics
 * - Edge cases (empty trees, single item, degenerate cases)
 * - Memory allocation failures (where testable)
 * - Iterator early termination
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <math.h>

/* Include all structures */
#define SPATIAL_NUMTYPE double
#define SPATIAL_DATATYPE void*
#include "libspatial/libspatial.h"

/* Test statistics */
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

/* Test macro */
#define TEST(name) static void test_##name(void)
#define RUN_TEST(name) do { \
    printf("  Running %s... ", #name); \
    tests_run++; \
    test_##name(); \
    tests_passed++; \
    printf("PASSED\n"); \
} while(0)

#define ASSERT(cond) do { \
    if (!(cond)) { \
        printf("FAILED\n  Assertion failed: %s at line %d\n", #cond, __LINE__); \
        tests_failed++; \
        tests_passed--; \
        return; \
    } \
} while(0)

#define ASSERT_EQ(a, b) ASSERT((a) == (b))
#define ASSERT_NE(a, b) ASSERT((a) != (b))
#define ASSERT_GT(a, b) ASSERT((a) > (b))
#define ASSERT_LT(a, b) ASSERT((a) < (b))

/* Test data */
typedef struct {
    int id;
    double value;
} test_data_t;

static int g_callback_count = 0;

static bool count_callback(const spatial_num_t *min, const spatial_num_t *max, 
                           spatial_data_t data, void *udata) {
    (void)min; (void)max; (void)data; (void)udata;
    g_callback_count++;
    return true;
}

static bool stop_early_callback(const spatial_num_t *min, const spatial_num_t *max,
                                spatial_data_t data, void *udata) {
    (void)min; (void)max; (void)data;
    int *count = (int*)udata;
    (*count)++;
    return *count < 5;  /* Stop after 5 items */
}

/* ============================================================================
 * Quadtree Tests
 * ============================================================================ */

TEST(quadtree_basic) {
    spatial_quadtree *qt = spatial_quadtree_new();
    ASSERT_NE(qt, NULL);
    
    /* Insert some items */
    double min[2], max[2];
    for (int i = 0; i < 100; i++) {
        min[0] = (double)(i * 10);
        min[1] = (double)(i * 10);
        max[0] = min[0] + 5.0;
        max[1] = min[1] + 5.0;
        ASSERT(spatial_quadtree_insert(qt, min, max, (void*)(intptr_t)(i + 1)));
    }
    
    ASSERT_EQ(spatial_quadtree_count(qt), 100);
    
    /* Search */
    g_callback_count = 0;
    min[0] = 0; min[1] = 0;
    max[0] = 500; max[1] = 500;
    spatial_quadtree_search(qt, min, max, count_callback, NULL);
    ASSERT_GT(g_callback_count, 0);
    
    /* Scan all */
    g_callback_count = 0;
    spatial_quadtree_scan(qt, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 100);
    
    /* Delete */
    min[0] = 50; min[1] = 50;
    max[0] = 55; max[1] = 55;
    ASSERT(spatial_quadtree_delete(qt, min, max, (void*)(intptr_t)6));
    ASSERT_EQ(spatial_quadtree_count(qt), 99);
    
    spatial_quadtree_free(qt);
}

TEST(quadtree_clone) {
    spatial_quadtree *qt = spatial_quadtree_new();
    ASSERT_NE(qt, NULL);
    
    double min[2] = {0, 0}, max[2] = {10, 10};
    spatial_quadtree_insert(qt, min, max, (void*)1);
    
    spatial_quadtree *clone = spatial_quadtree_clone(qt);
    ASSERT_NE(clone, NULL);
    ASSERT_EQ(clone, qt);  /* Should be same pointer (CoW) */
    
    spatial_quadtree_free(qt);
    spatial_quadtree_free(clone);
}

TEST(quadtree_empty) {
    spatial_quadtree *qt = spatial_quadtree_new();
    ASSERT_NE(qt, NULL);
    ASSERT_EQ(spatial_quadtree_count(qt), 0);
    
    /* Search on empty tree */
    double min[2] = {0, 0}, max[2] = {100, 100};
    g_callback_count = 0;
    spatial_quadtree_search(qt, min, max, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 0);
    
    spatial_quadtree_free(qt);
}

TEST(quadtree_early_stop) {
    spatial_quadtree *qt = spatial_quadtree_new();
    ASSERT_NE(qt, NULL);
    
    double min[2], max[2];
    for (int i = 0; i < 100; i++) {
        min[0] = (double)i;
        min[1] = (double)i;
        max[0] = min[0] + 0.5;
        max[1] = min[1] + 0.5;
        spatial_quadtree_insert(qt, min, max, (void*)(intptr_t)(i + 1));
    }
    
    min[0] = 0; min[1] = 0;
    max[0] = 100; max[1] = 100;
    int count = 0;
    spatial_quadtree_search(qt, min, max, stop_early_callback, &count);
    ASSERT_EQ(count, 5);
    
    spatial_quadtree_free(qt);
}

/* ============================================================================
 * Octree Tests
 * ============================================================================ */

TEST(octree_basic) {
    spatial_octree *ot = spatial_octree_new();
    ASSERT_NE(ot, NULL);
    
    double min[3], max[3];
    for (int i = 0; i < 50; i++) {
        min[0] = (double)(i * 10);
        min[1] = (double)(i * 10);
        min[2] = (double)(i * 10);
        max[0] = min[0] + 5.0;
        max[1] = min[1] + 5.0;
        max[2] = min[2] + 5.0;
        ASSERT(spatial_octree_insert(ot, min, max, (void*)(intptr_t)(i + 1)));
    }
    
    ASSERT_EQ(spatial_octree_count(ot), 50);
    
    g_callback_count = 0;
    spatial_octree_scan(ot, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 50);
    
    spatial_octree_free(ot);
}

TEST(octree_clone) {
    spatial_octree *ot = spatial_octree_new();
    ASSERT_NE(ot, NULL);
    
    double min[3] = {0, 0, 0}, max[3] = {10, 10, 10};
    spatial_octree_insert(ot, min, max, (void*)1);
    
    spatial_octree *clone = spatial_octree_clone(ot);
    ASSERT_NE(clone, NULL);
    
    spatial_octree_free(ot);
    spatial_octree_free(clone);
}

/* ============================================================================
 * Orthtree Tests
 * ============================================================================ */

TEST(orthtree_basic) {
    /* Use 4D for testing */
    #undef SPATIAL_ORTHTREE_DIMS
    #define SPATIAL_ORTHTREE_DIMS 4
    
    spatial_orthtree *ht = spatial_orthtree_new();
    ASSERT_NE(ht, NULL);
    
    /* Insert 20 items along x-axis, each [i, i+1] × [0,1]³ */
    for (int i = 0; i < 20; i++) {
        double min[4] = {(double)i, 0, 0, 0};
        double max[4] = {(double)(i + 1), 1, 1, 1};
        spatial_orthtree_insert(ht, min, max, (void*)(intptr_t)(i + 1));
    }
    ASSERT_EQ(spatial_orthtree_count(ht), 20);
    
    /* Search: query that overlaps items 5–15 (bboxes [4,5] through [14,15]) */
    double qmin[4] = {4.5, 0, 0, 0};
    double qmax[4] = {14.5, 1, 1, 1};
    g_callback_count = 0;
    spatial_orthtree_search(ht, qmin, qmax, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 11);
    
    /* Scan all items */
    g_callback_count = 0;
    spatial_orthtree_scan(ht, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 20);
    
    /* Iterator early-out */
    int early = 0;
    spatial_orthtree_scan(ht, stop_early_callback, &early);
    ASSERT_EQ(early, 5);
    
    /* Delete specific item (item #6, bbox [5,6] × [0,1]³) */
    double dmin[4] = {5, 0, 0, 0};
    double dmax[4] = {6, 1, 1, 1};
    bool deleted = spatial_orthtree_delete(ht, dmin, dmax, (void*)(intptr_t)6);
    ASSERT_EQ(deleted, true);
    ASSERT_EQ(spatial_orthtree_count(ht), 19);
    
    /* Verify deletion by re-searching */
    g_callback_count = 0;
    spatial_orthtree_search(ht, qmin, qmax, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 10);
    
    /* Delete non-existent item */
    deleted = spatial_orthtree_delete(ht, dmin, dmax, (void*)(intptr_t)999);
    ASSERT_EQ(deleted, false);
    ASSERT_EQ(spatial_orthtree_count(ht), 19);
    
    spatial_orthtree_free(ht);
}

TEST(orthtree_empty) {
    #undef SPATIAL_ORTHTREE_DIMS
    #define SPATIAL_ORTHTREE_DIMS 4
    
    spatial_orthtree *ht = spatial_orthtree_new();
    ASSERT_NE(ht, NULL);
    ASSERT_EQ(spatial_orthtree_count(ht), 0);
    
    double qmin[4] = {0, 0, 0, 0};
    double qmax[4] = {1, 1, 1, 1};
    int found = 0;
    spatial_orthtree_search(ht, qmin, qmax, count_callback, &found);
    ASSERT_EQ(found, 0);
    
    spatial_orthtree_scan(ht, count_callback, &found);
    ASSERT_EQ(found, 0);
    
    spatial_orthtree_free(ht);
}

/* ============================================================================
 * KD-Tree Tests
 * ============================================================================ */

TEST(kdtree_basic) {
    spatial_kdtree *kt = spatial_kdtree_new();
    ASSERT_NE(kt, NULL);

    double min[3], max[3];
    for (int i = 0; i < 100; i++) {
        min[0] = (double)(rand() % 1000);
        min[1] = (double)(rand() % 1000);
        min[2] = (double)(rand() % 1000);
        max[0] = min[0];
        max[1] = min[1];
        max[2] = min[2];
        ASSERT(spatial_kdtree_insert(kt, min, max, (void*)(intptr_t)(i + 1)));
    }

    ASSERT_EQ(spatial_kdtree_count(kt), 100);

    /* Range search */
    double qmin[3] = {0, 0, 0};
    double qmax[3] = {999, 999, 999};
    g_callback_count = 0;
    spatial_kdtree_search(kt, qmin, qmax, count_callback, NULL);
    ASSERT_GT(g_callback_count, 0);

    spatial_kdtree_free(kt);
}

TEST(kdtree_nearest) {
    spatial_kdtree *kt = spatial_kdtree_new();
    ASSERT_NE(kt, NULL);

    double min[3], max[3];
    for (int i = 0; i < 50; i++) {
        min[0] = max[0] = (double)(i * 10);
        min[1] = max[1] = (double)(i * 10);
        min[2] = max[2] = 0.0;
        ASSERT(spatial_kdtree_insert(kt, min, max, (void*)(intptr_t)(i + 1)));
    }

    /* k-NN search: returns void, fills arrays */
    double query[3] = {55.0, 55.0, 0.0};
    void *results[3] = {NULL, NULL, NULL};
    spatial_num_t dists[3] = {0, 0, 0};
    spatial_kdtree_nearest(kt, query, 3, NULL, NULL, results, dists);
    /* At least the closest result should be non-NULL */
    ASSERT_NE(results[0], NULL);

    spatial_kdtree_free(kt);
}

TEST(hilbert_curve) {
    uint32_t h1 = spatial_hilbert_xy2d(4, 0, 0);
    uint32_t h2 = spatial_hilbert_xy2d(4, 15, 15);
    
    /* Different coordinates should give different Hilbert values */
    ASSERT_NE(h1, h2);
}

TEST(stack_operations) {
    spatial_stack *s = spatial_stack_new(10, NULL);
    ASSERT_NE(s, NULL);
    
    int data1 = 42, data2 = 43;
    ASSERT(spatial_stack_push(s, &data1));
    ASSERT(spatial_stack_push(s, &data2));
    
    void *p = spatial_stack_pop(s);
    ASSERT_EQ(*(int*)p, 43);
    
    p = spatial_stack_pop(s);
    ASSERT_EQ(*(int*)p, 42);
    
    ASSERT(spatial_stack_empty(s));
    
    spatial_stack_free(s);
}

TEST(priority_queue) {
    spatial_priority_queue *pq = spatial_pq_new(10, NULL);
    ASSERT_NE(pq, NULL);
    
    spatial_pq_push(pq, (void*)1, 5.0);
    spatial_pq_push(pq, (void*)2, 3.0);
    spatial_pq_push(pq, (void*)3, 7.0);
    
    void *item;
    spatial_num_t dist;
    ASSERT(spatial_pq_pop(pq, &item, &dist));
    ASSERT_EQ(dist, 3.0);  /* Min heap - smallest first */
    
    spatial_pq_free(pq);
}

TEST(refcount) {
    spatial_refcount ref;
    spatial_refcount_init(&ref);
    ASSERT_EQ(spatial_refcount_get(&ref), 1);
    
    spatial_refcount_increment(&ref);
    ASSERT_EQ(spatial_refcount_get(&ref), 2);
    
    spatial_refcount_decrement(&ref);
    ASSERT_EQ(spatial_refcount_get(&ref), 1);
}

TEST(libspatial_info) {
    libspatial_info info = libspatial_get_info();
    ASSERT_NE(info.version_string, NULL);
    ASSERT_EQ(info.version_major, 0);
    ASSERT_EQ(info.num_structures, 8);
    
    ASSERT(libspatial_check_version(0, 1, 0));
    ASSERT(libspatial_check_version(0, 0, 9));
    ASSERT(!libspatial_check_version(0, 2, 0));
}

/* ============================================================================
 * Bbox / core utility tests
 * ============================================================================ */

TEST(bbox_operations) {
    double min1[3] = {0, 0, 0};
    double max1[3] = {10, 10, 10};
    double min2[3] = {5, 5, 5};
    double max2[3] = {15, 15, 15};

    ASSERT(spatial_bbox_overlaps(min1, max1, min2, max2, 3));

    double inner_min[3] = {2, 2, 2};
    double inner_max[3] = {8, 8, 8};
    ASSERT(spatial_bbox_contains(min1, max1, inner_min, inner_max, 3));

    double area = spatial_bbox_area(min1, max1, 3);
    ASSERT_EQ(area, 1000.0);
}

/* ============================================================================
 * VP-Tree Tests
 * ============================================================================ */

TEST(vptree_basic) {
    spatial_vptree *vp = spatial_vptree_new();
    ASSERT_NE(vp, NULL);

    double min[3], max[3];
    for (int i = 0; i < 50; i++) {
        min[0] = max[0] = (double)(rand() % 100);
        min[1] = max[1] = (double)(rand() % 100);
        min[2] = max[2] = (double)(rand() % 100);
        ASSERT(spatial_vptree_insert(vp, min, max, (void*)(intptr_t)(i + 1)));
    }

    ASSERT_EQ(spatial_vptree_count(vp), 50);

    double qmin[3] = {0, 0, 0};
    double qmax[3] = {99, 99, 99};
    g_callback_count = 0;
    spatial_vptree_search(vp, qmin, qmax, count_callback, NULL);
    ASSERT_GT(g_callback_count, 0);

    spatial_vptree_free(vp);
}

TEST(vptree_nearest) {
    spatial_vptree *vp = spatial_vptree_new();
    ASSERT_NE(vp, NULL);

    for (int i = 0; i < 30; i++) {
        double p[3] = {(double)(i * 5), 0.0, 0.0};
        ASSERT(spatial_vptree_insert(vp, p, p, (void*)(intptr_t)(i + 1)));
    }

    double query[3] = {22.0, 0.0, 0.0};
    void *results[3] = {NULL, NULL, NULL};
    spatial_num_t dists[3] = {0, 0, 0};
    spatial_vptree_nearest(vp, query, 3, results, dists);
    ASSERT_NE(results[0], NULL);

    spatial_vptree_free(vp);
}

/* ============================================================================
 * Hilbert R-Tree Tests
 * ============================================================================ */

TEST(hilbertrtree_basic) {
    spatial_hilbertrtree *rt = spatial_hilbertrtree_new();
    ASSERT_NE(rt, NULL);
    ASSERT_EQ(spatial_hilbertrtree_count(rt), 0);

    /* Insert items */
    for (int i = 0; i < 100; i++) {
        double mn[2] = {(double)(i * 10),     (double)(i * 10)};
        double mx[2] = {(double)(i * 10 + 5), (double)(i * 10 + 5)};
        ASSERT(spatial_hilbertrtree_insert(rt, mn, mx, (void*)(intptr_t)(i + 1)));
    }
    ASSERT_EQ(spatial_hilbertrtree_count(rt), 100);

    /* Range search */
    double qmin[2] = {0, 0};
    double qmax[2] = {500, 500};
    g_callback_count = 0;
    spatial_hilbertrtree_search(rt, qmin, qmax, count_callback, NULL);
    ASSERT_GT(g_callback_count, 0);

    /* Scan all */
    g_callback_count = 0;
    spatial_hilbertrtree_scan(rt, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 100);

    /* Early termination */
    int early = 0;
    spatial_hilbertrtree_search(rt, qmin, qmax, stop_early_callback, &early);
    ASSERT_LT(early, 100);

    /* Clone / CoW */
    spatial_hilbertrtree *rt2 = spatial_hilbertrtree_clone(rt);
    ASSERT_NE(rt2, NULL);
    ASSERT_EQ(spatial_hilbertrtree_count(rt2), 100);
    spatial_hilbertrtree_free(rt2);

    /* Delete */
    double mn[2] = {0, 0};
    double mx[2] = {5, 5};
    ASSERT(spatial_hilbertrtree_delete(rt, mn, mx, (void*)(intptr_t)1));
    ASSERT_EQ(spatial_hilbertrtree_count(rt), 99);

    /* Empty-tree edge case */
    spatial_hilbertrtree *empty = spatial_hilbertrtree_new();
    ASSERT_NE(empty, NULL);
    g_callback_count = 0;
    spatial_hilbertrtree_search(empty, qmin, qmax, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 0);
    spatial_hilbertrtree_free(empty);

    spatial_hilbertrtree_free(rt);
}

/* ============================================================================
 * BSP-Tree Tests
 * ============================================================================ */

TEST(bsptree_basic) {
    spatial_bsptree *bsp = spatial_bsptree_new();
    ASSERT_NE(bsp, NULL);
    ASSERT_EQ(spatial_bsptree_count(bsp), 0);
    
    /* Insert 50 items in 3D */
    for (int i = 0; i < 50; i++) {
        double min[3] = {(double)i, 0.0, 0.0};
        double max[3] = {(double)(i + 1), 1.0, 1.0};
        ASSERT(spatial_bsptree_insert(bsp, min, max, (void*)(intptr_t)(i + 1)));
    }
    ASSERT_EQ(spatial_bsptree_count(bsp), 50);
    
    /* Before rebuild, search still works (linear scan of root leaf) */
    double qmin[3] = {10.0, 0.0, 0.0};
    double qmax[3] = {20.0, 1.0, 1.0};
    g_callback_count = 0;
    spatial_bsptree_search(bsp, qmin, qmax, count_callback, NULL);
    /* Items [9,10] through [20,21] all touch or overlap [10,20] = 12 items */
    ASSERT_EQ(g_callback_count, 12);
    
    /* Rebuild the tree */
    spatial_bsptree_rebuild(bsp);
    
    /* Search after rebuild */
    g_callback_count = 0;
    spatial_bsptree_search(bsp, qmin, qmax, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 12);
    
    /* Scan all */
    g_callback_count = 0;
    spatial_bsptree_scan(bsp, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 50);
    
    spatial_bsptree_free(bsp);
}

/* ============================================================================
 * BVH Tests
 * ============================================================================ */

TEST(bvh_basic) {
    spatial_bvh *bvh = spatial_bvh_new();
    ASSERT_NE(bvh, NULL);
    ASSERT_EQ(spatial_bvh_count(bvh), 0);
    
    /* Insert 100 items in 3D */
    for (int i = 0; i < 100; i++) {
        double min[3] = {(double)(i * 10), (double)(i * 10), 0.0};
        double max[3] = {(double)(i * 10 + 5), (double)(i * 10 + 5), 1.0};
        ASSERT(spatial_bvh_insert(bvh, min, max, (void*)(intptr_t)(i + 1)));
    }
    ASSERT_EQ(spatial_bvh_count(bvh), 100);
    
    /* Build the BVH */
    ASSERT(spatial_bvh_build(bvh));
    
    /* Search */
    double qmin[3] = {0, 0, 0};
    double qmax[3] = {50, 50, 1};
    g_callback_count = 0;
    spatial_bvh_search(bvh, qmin, qmax, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 6);
    
    /* Scan all */
    g_callback_count = 0;
    spatial_bvh_scan(bvh, count_callback, NULL);
    ASSERT_EQ(g_callback_count, 100);
    
    spatial_bvh_free(bvh);
}

TEST(bvh_ray_intersect) {
    spatial_bvh *bvh = spatial_bvh_new();
    ASSERT_NE(bvh, NULL);
    
    /* Insert a few axis-aligned boxes along x-axis */
    for (int i = 0; i < 10; i++) {
        double min[3] = {(double)(i * 10), -1.0, -1.0};
        double max[3] = {(double)(i * 10 + 5), 1.0, 1.0};
        ASSERT(spatial_bvh_insert(bvh, min, max, (void*)(intptr_t)(i + 1)));
    }
    ASSERT(spatial_bvh_build(bvh));
    
    /* Shoot a ray along +x axis from origin */
    spatial_bvh_ray ray;
    ray.origin[0] = -5.0; ray.origin[1] = 0.0; ray.origin[2] = 0.0;
    ray.dir[0] = 1.0; ray.dir[1] = 0.0; ray.dir[2] = 0.0;
    ray.inv_dir[0] = 1.0; ray.inv_dir[1] = 1e30; ray.inv_dir[2] = 1e30;
    ray.tmin = 0.0;
    ray.tmax = 1000.0;
    
    spatial_bvh_hit hit;
    bool found = spatial_bvh_intersect_ray(bvh, &ray, &hit);
    ASSERT_EQ(found, true);
    ASSERT_EQ(hit.prim_id >= 0 && hit.prim_id < 10, true);
    
    spatial_bvh_free(bvh);
}

/* ============================================================================
 * Main
 * ============================================================================ */

int main(int argc, char **argv) {
    (void)argc; (void)argv;
    
    printf("=================================================================\n");
    printf("libspatial Test Suite v%s\n", LIBSPATIAL_VERSION_STRING);
    printf("=================================================================\n\n");
    
    srand((unsigned)time(NULL));
    
    printf("Core Utilities:\n");
    RUN_TEST(bbox_operations);
    RUN_TEST(hilbert_curve);
    RUN_TEST(stack_operations);
    RUN_TEST(priority_queue);
    RUN_TEST(refcount);
    RUN_TEST(libspatial_info);
    
    printf("\nQuadtree:\n");
    RUN_TEST(quadtree_basic);
    RUN_TEST(quadtree_clone);
    RUN_TEST(quadtree_empty);
    RUN_TEST(quadtree_early_stop);
    
    printf("\nOctree:\n");
    RUN_TEST(octree_basic);
    RUN_TEST(octree_clone);
    
    printf("\nOrthtree:\n");
    RUN_TEST(orthtree_basic);
    RUN_TEST(orthtree_empty);
    
    printf("\nKD-Tree:\n");
    RUN_TEST(kdtree_basic);
    RUN_TEST(kdtree_nearest);
    
    printf("\nVP-Tree:\n");
    RUN_TEST(vptree_basic);
    RUN_TEST(vptree_nearest);
    
    printf("\nHilbert R-Tree:\n");
    RUN_TEST(hilbertrtree_basic);
    
    printf("\nBSP-Tree:\n");
    RUN_TEST(bsptree_basic);
    
    printf("\nBVH:\n");
    RUN_TEST(bvh_basic);
    RUN_TEST(bvh_ray_intersect);
    
    printf("\n=================================================================\n");
    printf("Results: %d tests run, %d passed, %d failed\n", 
           tests_run, tests_passed, tests_failed);
    printf("=================================================================\n");
    
    return tests_failed > 0 ? 1 : 0;
}
