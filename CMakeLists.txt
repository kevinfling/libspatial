# -----------------------------------------------------------------------------
# CMakeLists.txt
# Build configuration for libspatial
# Version: 0.1.0
# Author:  Kevin Fling
# License: MIT
# Copyright (c) 2026 Kevin Fling
#
# libspatial — Header-only C11 library providing production-grade spatial
# indexing data structures.
# -----------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.14)
project(libspatial VERSION 0.1.0 LANGUAGES C)

# C11 is required
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Optional memento integration (advanced arena/pool allocator)
option(LIBSPATIAL_USE_MEMENTO "Fetch and enable memento integration" ON)

if(LIBSPATIAL_USE_MEMENTO)
    include(FetchContent)
    FetchContent_Declare(
        memento
        GIT_REPOSITORY https://github.com/StraylightRunMOO/memento.git
        GIT_TAG        master
    )
    FetchContent_GetProperties(memento)
    if(NOT memento_POPULATED)
        FetchContent_Populate(memento)
    endif()
endif()

# Interface library
add_library(libspatial INTERFACE)
add_library(libspatial::libspatial ALIAS libspatial)

target_include_directories(libspatial INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Compiler flags for consumers (strict)
target_compile_options(libspatial INTERFACE
    $<$<C_COMPILER_ID:GNU,Clang,AppleClang>:
        -std=c11
        -pedantic
        -Wall
        -Wextra
        -Werror
        -Wstrict-prototypes
        -Wmissing-prototypes
        -Wshadow
        -Wcast-qual
        -Wpointer-arith
        -Wwrite-strings
        -Wswitch-enum
        -Winit-self
        -Wmissing-field-initializers
        -Wredundant-decls
        -Wfloat-equal
        -Wundef
        -Wvla
    >
    $<$<C_COMPILER_ID:MSVC>:
        /W4
        /WX
        /permissive-
        /wd4996
    >
)

# Options
option(BUILD_TESTING "Build test suite" ON)
option(BUILD_BENCHMARKS "Build benchmarks" ON)
option(BUILD_EXAMPLES "Build examples" ON)

# Testing
if(BUILD_TESTING)
    enable_testing()
    add_executable(test_libspatial tests/test_libspatial.c)
    target_link_libraries(test_libspatial PRIVATE libspatial)
    target_compile_definitions(test_libspatial PRIVATE SPATIAL_TESTING)
    
    # Additional sanitizer flags for tests
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(test_libspatial PRIVATE
            -fsanitize=address,undefined
            -fno-omit-frame-pointer
            -g
            -Wno-static-in-inline
        )
        target_link_options(test_libspatial PRIVATE
            -fsanitize=address,undefined
        )
    endif()
    
    add_test(NAME libspatial_test COMMAND test_libspatial)
endif()

# Benchmarks
if(BUILD_BENCHMARKS)
    add_executable(bench_libspatial benchmarks/bench_libspatial.c)
    target_link_libraries(bench_libspatial PRIVATE libspatial)
    if(LIBSPATIAL_USE_MEMENTO)
        target_include_directories(bench_libspatial PRIVATE ${memento_SOURCE_DIR}/include)
        target_compile_definitions(bench_libspatial PRIVATE SPATIAL_BENCHMARKING MEMENTO_IMPLEMENTATION)
    else()
        target_compile_definitions(bench_libspatial PRIVATE SPATIAL_BENCHMARKING)
    endif()
    
    # Optimize for benchmarks
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
        target_compile_options(bench_libspatial PRIVATE
            -O3
            -march=native
            -ffast-math
            -funroll-loops
            -Wno-static-in-inline
        )
    elseif(MSVC)
        target_compile_options(bench_libspatial PRIVATE /O2 /arch:AVX2)
    endif()
endif()

# Examples
if(BUILD_EXAMPLES)
    add_executable(arena_example examples/arena_example.c)
    target_link_libraries(arena_example PRIVATE libspatial m)

    add_executable(city_search examples/city_search.c)
    target_link_libraries(city_search PRIVATE libspatial m)

    if(LIBSPATIAL_USE_MEMENTO)
        add_executable(memento_arena_example examples/memento_arena_example.c)
        target_link_libraries(memento_arena_example PRIVATE libspatial m)
        target_include_directories(memento_arena_example PRIVATE ${memento_SOURCE_DIR}/include)
        target_compile_definitions(memento_arena_example PRIVATE MEMENTO_IMPLEMENTATION)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(memento_arena_example PRIVATE -Wno-static-in-inline -Wno-cast-align)
        endif()

        add_executable(memento_pool_example examples/memento_pool_example.c)
        target_link_libraries(memento_pool_example PRIVATE libspatial m)
        target_include_directories(memento_pool_example PRIVATE ${memento_SOURCE_DIR}/include)
        target_compile_definitions(memento_pool_example PRIVATE MEMENTO_IMPLEMENTATION)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
            target_compile_options(memento_pool_example PRIVATE -Wno-static-in-inline -Wno-cast-align)
        endif()
    endif()
endif()

# Installation
include(GNUInstallDirs)

install(TARGETS libspatial
    EXPORT libspatial-targets
)

install(EXPORT libspatial-targets
    FILE libspatial-targets.cmake
    NAMESPACE libspatial::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libspatial
)

install(DIRECTORY include/libspatial
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Package config
include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/libspatial-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/libspatial-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libspatial
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/libspatial-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/libspatial-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/libspatial-config-version.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/libspatial
)

# Summary
message(STATUS "")
message(STATUS "libspatial ${PROJECT_VERSION} configuration:")
message(STATUS "  Use memento:      ${LIBSPATIAL_USE_MEMENTO}")
message(STATUS "  Build testing:    ${BUILD_TESTING}")
message(STATUS "  Build benchmarks: ${BUILD_BENCHMARKS}")
message(STATUS "  Build examples:   ${BUILD_EXAMPLES}")
message(STATUS "  Install prefix:   ${CMAKE_INSTALL_PREFIX}")
message(STATUS "")