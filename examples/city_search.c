/**
 * @file   city_search.c
 * @brief  Example: City search using Hilbert R-Tree
 * @version 0.1.0
 * @author Kevin Fling
 * @license MIT
 * @copyright Copyright (c) 2026 Kevin Fling
 *
 * libspatial — Header-only C11 library providing production-grade spatial
 * indexing data structures.
 *
 * This example demonstrates using libspatial for geographic searches,
 * similar to the city-search example in tidwall/rtree.c
 * Usage:
 * ./city_search "New York"
 * ./city_search -near 40.7128 -74.0060 100  # Cities within 100km of NYC
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "libspatial/libspatial.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Sample city database */
typedef struct {
    const char *name;
    double lat;  /* Latitude */
    double lon;  /* Longitude */
    int population;
} city_t;

static city_t cities[] = {
    {"New York", 40.7128, -74.0060, 8400000},
    {"Los Angeles", 34.0522, -118.2437, 3900000},
    {"Chicago", 41.8781, -87.6298, 2700000},
    {"Houston", 29.7604, -95.3698, 2300000},
    {"Phoenix", 33.4484, -112.0740, 1600000},
    {"Philadelphia", 39.9526, -75.1652, 1500000},
    {"San Antonio", 29.4241, -98.4936, 1400000},
    {"San Diego", 32.7157, -117.1611, 1400000},
    {"Dallas", 32.7767, -96.7970, 1300000},
    {"San Jose", 37.3382, -121.8863, 1000000},
    {"Austin", 30.2672, -97.7431, 950000},
    {"Jacksonville", 30.3322, -81.6557, 900000},
    {"Fort Worth", 32.7555, -97.3308, 890000},
    {"Columbus", 39.9612, -82.9988, 880000},
    {"Charlotte", 35.2271, -80.8431, 870000},
    {"San Francisco", 37.7749, -122.4194, 870000},
    {"Indianapolis", 39.7684, -86.1581, 860000},
    {"Seattle", 47.6062, -122.3321, 750000},
    {"Denver", 39.7392, -104.9903, 715000},
    {"Washington", 38.9072, -77.0369, 700000},
    {"Boston", 42.3601, -71.0589, 690000},
    {"El Paso", 31.7619, -106.4850, 680000},
    {"Nashville", 36.1627, -86.7816, 670000},
    {"Detroit", 42.3314, -83.0458, 670000},
    {"Oklahoma City", 35.4676, -97.5164, 650000},
    {"Portland", 45.5152, -122.6784, 650000},
    {"Las Vegas", 36.1699, -115.1398, 640000},
    {"Louisville", 38.2527, -85.7585, 620000},
    {"Baltimore", 39.2904, -76.6122, 610000},
    {"Milwaukee", 43.0389, -87.9065, 590000},
    {NULL, 0, 0, 0}
};

/* Convert lat/lon to approximate meters (equirectangular projection) */
static void latlon_to_meters(double lat, double lon, double *x, double *y) {
    const double R = 6371000.0;  /* Earth radius in meters */
    const double lat0 = 0.0;     /* Reference latitude */
    
    *x = R * lon * M_PI / 180.0 * cos(lat0 * M_PI / 180.0);
    *y = R * lat * M_PI / 180.0;
}

/* Haversine distance in km */
static double haversine_distance(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371.0;  /* Earth radius in km */
    
    double dlat = (lat2 - lat1) * M_PI / 180.0;
    double dlon = (lon2 - lon1) * M_PI / 180.0;
    
    double a = sin(dlat / 2) * sin(dlat / 2) +
               cos(lat1 * M_PI / 180.0) * cos(lat2 * M_PI / 180.0) *
               sin(dlon / 2) * sin(dlon / 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    
    return R * c;
}

/* Search context */
typedef struct {
    const char *query;
    double lat, lon;
    double radius_km;
    int found;
    int mode;  /* 0 = name search, 1 = radius search */
} search_ctx_t;

/* Portable case-insensitive substring search */
static const char *str_case_str(const char *haystack, const char *needle) {
    if (!needle[0]) return haystack;
    for (const char *h = haystack; *h; h++) {
        if (tolower((unsigned char)*h) == tolower((unsigned char)*needle)) {
            const char *h2 = h + 1;
            const char *n2 = needle + 1;
            while (*n2 && tolower((unsigned char)*h2) == tolower((unsigned char)*n2)) {
                h2++; n2++;
            }
            if (!*n2) return h;
        }
    }
    return NULL;
}

/* Callback for name search */
static bool name_search_callback(const spatial_num_t *min, const spatial_num_t *max,
                                  spatial_data_t data, void *udata) {
    (void)min; (void)max;
    search_ctx_t *ctx = (search_ctx_t*)udata;
    city_t *city = (city_t*)data;
    
    if (str_case_str(city->name, ctx->query)) {
        printf("  %s (pop: %d, lat: %.4f, lon: %.4f)\n",
               city->name, city->population, city->lat, city->lon);
        ctx->found++;
    }
    
    return true;
}

/* Callback for radius search */
static bool radius_search_callback(const spatial_num_t *min, const spatial_num_t *max,
                                    spatial_data_t data, void *udata) {
    (void)min; (void)max;
    search_ctx_t *ctx = (search_ctx_t*)udata;
    city_t *city = (city_t*)data;
    
    double dist = haversine_distance(ctx->lat, ctx->lon, city->lat, city->lon);
    
    if (dist <= ctx->radius_km) {
        printf("  %s (%.1f km away, pop: %d)\n", city->name, dist, city->population);
        ctx->found++;
    }
    
    return true;
}

static void print_usage(const char *prog) {
    printf("Usage: %s [options] <query>\n", prog);
    printf("\n");
    printf("Options:\n");
    printf("  -near <lat> <lon> <radius_km>  Find cities within radius\n");
    printf("  -h, --help                      Show this help\n");
    printf("\n");
    printf("Examples:\n");
    printf("  %s \"San\"              # Cities containing 'San'\n", prog);
    printf("  %s -near 40.7 -74.0 100  # Cities within 100km of NYC\n", prog);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }
    
    /* Parse arguments */
    search_ctx_t ctx = {0};
    
    if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
        print_usage(argv[0]);
        return 0;
    }
    
    if (strcmp(argv[1], "-near") == 0) {
        if (argc < 6) {
            fprintf(stderr, "Error: -near requires lat, lon, and radius\n");
            return 1;
        }
        ctx.mode = 1;
        ctx.lat = atof(argv[2]);
        ctx.lon = atof(argv[3]);
        ctx.radius_km = atof(argv[4]);
    } else {
        ctx.mode = 0;
        ctx.query = argv[1];
    }
    
    /* Build spatial index */
    printf("Building spatial index...\n");
    
    spatial_hilbertrtree *rt = spatial_hilbertrtree_new();
    if (!rt) {
        fprintf(stderr, "Error: Failed to create spatial index\n");
        return 1;
    }
    
    for (int i = 0; cities[i].name; i++) {
        double x, y;
        latlon_to_meters(cities[i].lat, cities[i].lon, &x, &y);
        
        /* Small bounding box around each city point */
        double min[2] = {x - 1000, y - 1000};  /* 1km buffer */
        double max[2] = {x + 1000, y + 1000};
        
        spatial_hilbertrtree_insert(rt, min, max, &cities[i]);
    }
    
    printf("Indexed %d cities\n\n", spatial_hilbertrtree_count(rt));
    
    /* Perform search */
    if (ctx.mode == 0) {
        /* Name search */
        printf("Searching for cities containing '%s':\n", ctx.query);
        
        /* Search entire world bounds */
        double qmin[2] = {-1e9, -1e9};
        double qmax[2] = {1e9, 1e9};
        
        spatial_hilbertrtree_search(rt, qmin, qmax, name_search_callback, &ctx);
        
        if (ctx.found == 0) {
            printf("  No cities found\n");
        } else {
            printf("\nFound %d cities\n", ctx.found);
        }
    } else {
        /* Radius search */
        printf("Searching for cities within %.1f km of (%.4f, %.4f):\n",
               ctx.radius_km, ctx.lat, ctx.lon);
        
        /* Convert search center to meters and create query box */
        double cx, cy;
        latlon_to_meters(ctx.lat, ctx.lon, &cx, &cy);
        
        /* Query box (slightly larger than radius to catch edge cases) */
        double r_m = ctx.radius_km * 1000.0 * 1.5;
        double qmin[2] = {cx - r_m, cy - r_m};
        double qmax[2] = {cx + r_m, cy + r_m};
        
        spatial_hilbertrtree_search(rt, qmin, qmax, radius_search_callback, &ctx);
        
        if (ctx.found == 0) {
            printf("  No cities found within radius\n");
        } else {
            printf("\nFound %d cities within %.1f km\n", ctx.found, ctx.radius_km);
        }
    }
    
    spatial_hilbertrtree_free(rt);
    
    return 0;
}
