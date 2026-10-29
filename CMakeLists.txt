cmake_minimum_required(VERSION 3.24)
project(pavel LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(cmake/Dependencies.cmake)

add_subdirectory(src)
add_subdirectory(app)

enable_testing()
add_subdirectory(tests)
