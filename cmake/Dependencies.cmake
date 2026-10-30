include(FetchContent)

# libigl comes first: its Eigen recipe defines Eigen3::Eigen, which geometry-central then reuses.
FetchContent_Declare(libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG v2.6.4
    GIT_SHALLOW TRUE)

FetchContent_Declare(geometry-central
    GIT_REPOSITORY https://github.com/nmwsharp/geometry-central.git
    GIT_TAG v1.1.0
    GIT_SHALLOW TRUE
    SYSTEM) # its headers trip our -Wextra

# Release tarball instead of git: it already contains Geogram's submodules.
FetchContent_Declare(geogram
    URL https://github.com/BrunoLevy/geogram/releases/download/v1.10.1/geogram_1.10.1.tar.gz
    URL_HASH SHA256=30033340a0dfd86f9e15a8ffdd7c1e6d49b967d3f26fb96cf3dfd1a502d8eba2)
set(GEOGRAM_LIB_ONLY ON)
set(GEOGRAM_WITH_GRAPHICS OFF)
set(GEOGRAM_WITH_LUA OFF)

FetchContent_Declare(CLI11
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG v2.7.2
    GIT_SHALLOW TRUE)

FetchContent_Declare(Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.16.1
    GIT_SHALLOW TRUE)

FetchContent_MakeAvailable(libigl geometry-central geogram CLI11 Catch2)

option(PAVEL_VIEWER "Build the polyscope viewer for manual placement (pavel --place)" OFF)
if(PAVEL_VIEWER)
    FetchContent_Declare(polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG v2.6.1
        GIT_SHALLOW TRUE
        SYSTEM)
    FetchContent_MakeAvailable(polyscope)
endif()
