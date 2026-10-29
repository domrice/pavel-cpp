add_executable(pavel_tests test_recovery_function.cpp test_voxelize.cpp)
target_link_libraries(pavel_tests PRIVATE pavel::core Catch2::Catch2WithMain)

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(pavel_tests)
