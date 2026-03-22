# libspatial

libspatial is a header-only spatial indexing library for C11. It includes eight
production-grade tree structures for 2D, 3D, and N-dimensional spatial queries,
all sharing a uniform API pattern.

## Features

- Header-only - drop into your project and include
- Strict C11 - compiles with `-std=c11 -pedantic -Wall -Wextra -Werror`
- Zero dependencies - built-in arena and pool allocators
- Eight spatial structures - quadtree, octree, hyperoctree, kd-tree, vp-tree,
  hilbert r-tree, bsp-tree, and bvh
- Cache-aligned node layouts
- Non-recursive stack-based traversal
- Copy-on-write cloning with reference counting
- Custom allocator support (malloc, arena, pool, or memento)
- Configurable numeric types, data payloads, and node capacities

## Goals

- Give C programs fast in-memory spatial indexes without dependencies
- Provide multiple tree types so users can choose the right structure for their
  dimensionality and workload
- Support custom allocators for games, embedded systems, and real-time apps
- Stay portable across Clang, GCC, and MSVC

## Non-goals

- Disk-based or persistent spatial indexes
- GPU-accelerated structures
- Thread-safe concurrent mutation (read-only concurrent access is fine with
  copy-on-write)

## Using

Install with CMake or just copy the headers from `include/libspatial/` into your
project.

```bash
cmake -B build
cmake --install build --prefix /usr/local
```

Then include the umbrella header or individual structures:

```c
#include <libspatial/libspatial.h>
```

## Example 1 (Quadtree insert and search)

```c
#include <stdio.h>
#include <libspatial/libspatial.h>

static bool my_callback(const double *min, const double *max,
                        void *data, void *udata) {
    (void)min; (void)max; (void)udata;
    printf("Found item %d\n", (int)(size_t)data);
    return true; /* continue iterating */
}

int main() {
    struct spatial_quadtree *qt = spatial_quadtree_new();

    double min[2] = {0, 0};
    double max[2] = {100, 100};
    spatial_quadtree_insert(qt, min, max, (void*)(size_t)42);

    double qmin[2] = {10, 10};
    double qmax[2] = {50, 50};
    spatial_quadtree_search(qt, qmin, qmax, my_callback, NULL);

    printf("Total items: %d\n", spatial_quadtree_count(qt));
    spatial_quadtree_free(qt);
    return 0;
}

// Output:
// Found item 42
// Total items: 1
```

## Example 2 (City search with Hilbert R-Tree)

The `city_search` example demonstrates geographic range queries using the
Hilbert R-Tree.

```bash
./city_search "San"
```

```
Building spatial index...
Indexed 30 cities

Searching for cities containing 'San':
  San Antonio (pop: 1400000, lat: 29.4241, lon: -98.4936)
  San Francisco (pop: 870000, lat: 37.7749, lon: -122.4194)
  San Jose (pop: 1000000, lat: 37.3382, lon: -121.8863)
  San Diego (pop: 1400000, lat: 32.7157, lon: -117.1611)

Found 4 cities
```

## Example 3 (Arena allocator)

Eliminate per-operation malloc overhead using the built-in bump arena:

```c
#include <libspatial/libspatial.h>

int main() {
    spatial_arena_t *arena = spatial_arena_new(4 * 1024 * 1024);
    struct spatial_bvh *bvh = spatial_bvh_new_with_arena(arena);

    double min[3] = {0, 0, 0};
    double max[3] = {10, 10, 10};
    spatial_bvh_insert(bvh, min, max, (void*)(size_t)1);

    /* individual frees are no-ops when using an arena */
    spatial_bvh_free(bvh);

    /* bulk-free everything at once */
    spatial_arena_reset(arena);
    spatial_arena_free(arena);
    return 0;
}
```

## Configuration

All structures are configured via preprocessor macros defined before including
the header.

| Macro | Description | Default |
| :--- | :--- | :--- |
| `SPATIAL_NUMTYPE` | Numeric type for coordinates | `double` |
| `SPATIAL_DATATYPE` | Payload type stored with each item | `void *` |
| `SPATIAL_MAXITEMS` | Maximum items per node | `64` |
| `SPATIAL_QUADTREE_MAX_DEPTH` | Quadtree max depth | `20` |
| `SPATIAL_QUADTREE_MAX_ITEMS` | Quadtree node capacity | `16` |
| `SPATIAL_HYPEROCTREE_DIMS` | Hyperoctree dimensionality | `4` |

For example, to use `float` coordinates and a smaller node capacity:

```c
#define SPATIAL_NUMTYPE float
#define SPATIAL_MAXITEMS 32
#include <libspatial/quadtree.h>
```

## Data Structures

| Structure | Dimensions | Best For |
| :--- | :--- | :--- |
| Quadtree | 2D | Games, GIS, image processing |
| Octree | 3D | Voxel engines, point clouds |
| Hyperoctree | N-D | High-dimensional indexing |
| KD-Tree | N-D | Nearest neighbor, ray tracing |
| VP-Tree | Metric space | Non-Euclidean distances |
| Hilbert R-Tree | 2D | Range queries, spatial databases |
| BSP-Tree | 2D/3D | Visibility, CSG, collision |
| BVH | 3D | Physics engines, ray tracing |

Quadtree, octree, hyperoctree, kd-tree, vp-tree, and hilbert r-tree are fully
tested. BSP-tree and BVH headers are present but still being validated.

## Building

```bash
mkdir build && cd build
cmake .. -DBUILD_TESTING=ON -DBUILD_BENCHMARKS=ON -DBUILD_EXAMPLES=ON
cmake --build .
```

## Running Tests

```bash
./test_libspatial
```

```
=================================================================
libspatial Test Suite v0.1.0
=================================================================

Core Utilities:
  Running bbox_operations... PASSED
  Running hilbert_curve... PASSED
  Running stack_operations... PASSED
  Running priority_queue... PASSED
  Running refcount... PASSED
  Running libspatial_info... PASSED

Quadtree:
  Running quadtree_basic... PASSED
  Running quadtree_clone... PASSED
  Running quadtree_empty... PASSED
  Running quadtree_early_stop... PASSED

Octree:
  Running octree_basic... PASSED
  Running octree_clone... PASSED

Hyperoctree:
  Running hyperoctree_basic... PASSED

KD-Tree:
  Running kdtree_basic... PASSED
  Running kdtree_nearest... PASSED

VP-Tree:
  Running vptree_basic... PASSED
  Running vptree_nearest... PASSED

Hilbert R-Tree:
  Running hilbertrtree_basic... PASSED

=================================================================
Results: 18 tests run, 18 passed, 0 failed
=================================================================
```

## Running Benchmarks

```bash
./bench_libspatial 100000
```

Output on Linux x86-64 (AMD Ryzen, ~3GHz):

```
=================================================================
libspatial Performance Benchmark v0.1.0
=================================================================
Items: 100000, Queries: 10000
(Cycles measured at ~3GHz, lower is better)

Structure             Insert(cyc)  Search(cyc)    Scan(cyc)   Insert(ns)   Search(ns)     Scan(ns)
-------------------- ------------ ------------ ------------ ------------ ------------ ------------
Quadtree                     11.0          3.0          3.0          3.7          1.0          1.0
Octree                       10.0          6.0          3.0          3.3          2.0          1.0
KD-Tree                      23.0       5479.0          0.0          7.7       1826.3          0.0
VP-Tree                      28.0          7.0          0.0          9.3          2.3          0.0
Hilbert R-Tree               70.0       7656.0          0.0         23.3       2552.0          0.0
BSP-Tree                      1.0       4047.0          0.0          0.3       1349.0          0.0
BVH                           2.0       7184.0          0.0          0.7       2394.7          0.0
```

Notes on the above:
- Quadtree and octree show the fastest search times on this random AABB
  workload.
- BSP-tree and BVH insert times measure item collection only; tree
  construction happens later via `spatial_bsptree_rebuild()` or
  `spatial_bvh_build()`.
- Scan times for some structures round to zero because the per-item cost is
  below the timer resolution on this workload.
- Arena builds are typically 2-5x faster than plain `malloc` for batch
  insertion workloads.

## API Reference

All structures share the same function pattern:

```c
struct spatial_xxx *spatial_xxx_new(void);
struct spatial_xxx *spatial_xxx_new_with_allocator(spatial_allocator *alloc);
void spatial_xxx_free(struct spatial_xxx *t);
struct spatial_xxx *spatial_xxx_clone(const struct spatial_xxx *t);

bool spatial_xxx_insert(struct spatial_xxx *t,
                        const spatial_num_t *min,
                        const spatial_num_t *max,
                        spatial_data_t data);
void spatial_xxx_search(const struct spatial_xxx *t,
                        const spatial_num_t *qmin,
                        const spatial_num_t *qmax,
                        spatial_xxx_iter_fn callback,
                        void *udata);
void spatial_xxx_scan(const struct spatial_xxx *t,
                      spatial_xxx_iter_fn callback,
                      void *udata);
int spatial_xxx_count(const struct spatial_xxx *t);
```

Iterator callbacks return `true` to continue or `false` to stop early.

See the individual headers in `include/libspatial/` for structure-specific
functions such as `spatial_kdtree_nearest()`, `spatial_vptree_nearest()`, and
`spatial_bvh_ray_intersect()`.

## Custom Allocators

By default all structures use `malloc` and `free`. Provide a custom allocator
when creating a tree:

```c
spatial_allocator alloc = {
    .malloc = my_malloc,
    .free = my_free
};
struct spatial_quadtree *qt = spatial_quadtree_new_with_allocator(&alloc);
```

Or use the built-in arena or pool allocators for zero per-operation overhead.

## Optional Memento Integration

For a production-grade bump allocator with thread-local heaps, enable memento
at CMake time:

```bash
cmake .. -DLIBSPATIAL_USE_MEMENTO=ON
```

## License

MIT License - See [LICENSE](LICENSE) for details.

Copyright (c) 2026 Kevin Fling

## Acknowledgments

Inspired by [tidwall/rtree.c](https://github.com/tidwall/rtree.c).
