add_executable(pavel main.cpp)
target_link_libraries(pavel PRIVATE pavel::core CLI11::CLI11)
