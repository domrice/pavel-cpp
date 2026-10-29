add_executable(pavel_tests test_fast_marching.cpp test_recovery_function.cpp test_voxelize.cpp)
target_link_libraries(pavel_tests PRIVATE pavel::core Catch2::Catch2WithMain)
target_compile_definitions(pavel_tests PRIVATE PAVEL_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/data")

list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
include(Catch)
catch_discover_tests(pavel_tests)
