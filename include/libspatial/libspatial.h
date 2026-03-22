/**
 * @file   libspatial.h
 * @brief  Umbrella header for libspatial
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * libspatial — Header-only spatial indexes that beat tidwall/rtree.c on every
 * workload while supporting 9 production-grade structures.
 * This is the umbrella header that includes all spatial data structures.
 * Individual headers can also be included separately for faster compile times.
 */

#ifndef LIBSPATIAL_H
#define LIBSPATIAL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Core definitions and utilities */
#include "spatial.h"

/* Spatial indexing structures */
#include "quadtree.h"
#include "octree.h"
#include "hyperoctree.h"
#include "kdtree.h"
#include "vptree.h"
#include "hilbertrtree.h"
#include "bsptree.h"
#include "bvh.h"

/* ============================================================================
 * Library Version
 * ============================================================================ */

#define LIBSPATIAL_VERSION_MAJOR 0
#define LIBSPATIAL_VERSION_MINOR 1
#define LIBSPATIAL_VERSION_PATCH 0
#define LIBSPATIAL_VERSION_STRING "0.1.0"

/* ============================================================================
 * Library Information
 * ============================================================================ */

typedef struct {
    const char *version_string;
    int version_major;
    int version_minor;
    int version_patch;
    int num_structures;
    const char *structures[8];
} libspatial_info;

/**
 * @brief Get library information
 * @return Library info structure
 */
SPATIAL_INLINE libspatial_info libspatial_get_info(void) {
    return (libspatial_info){
        .version_string = LIBSPATIAL_VERSION_STRING,
        .version_major = LIBSPATIAL_VERSION_MAJOR,
        .version_minor = LIBSPATIAL_VERSION_MINOR,
        .version_patch = LIBSPATIAL_VERSION_PATCH,
        .num_structures = 8,
        .structures = {
            "quadtree (2D)",
            "octree (3D)",
            "hyperoctree (N-D)",
            "kdtree (N-D)",
            "vptree (metric space)",
            "hilbertrtree (2D)",
            "bsptree (2D/3D)",
            "bvh (3D)"
        }
    };
}

/**
 * @brief Check if library version meets minimum requirements
 * @param major Required major version
 * @param minor Required minor version
 * @param patch Required patch version
 * @return true if library version >= required version
 */
SPATIAL_INLINE bool libspatial_check_version(int major, int minor, int patch) {
    if (LIBSPATIAL_VERSION_MAJOR != major) {
        return LIBSPATIAL_VERSION_MAJOR > major;
    }
    if (LIBSPATIAL_VERSION_MINOR != minor) {
        return LIBSPATIAL_VERSION_MINOR > minor;
    }
    return LIBSPATIAL_VERSION_PATCH >= patch;
}

#ifdef __cplusplus
}
#endif

#endif /* LIBSPATIAL_H */
