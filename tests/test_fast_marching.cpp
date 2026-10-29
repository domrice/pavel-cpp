#include "fast_marching.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <fstream>
#include <iterator>
#include <string>

// Data of a C-order .npy (v1 header) as an Eigen array. The fixtures are 32³ indexed [i, j, k]
// with k fastest, which is our x-fastest storage with x = k.
template <class T>
Eigen::Array<T, Eigen::Dynamic, 1> load_npy(const std::string& name)
{
    std::ifstream f(std::string(PAVEL_TEST_DATA) + "/fmm/" + name + ".npy", std::ios::binary);
    const std::string bytes(std::istreambuf_iterator<char>(f), {});
    const auto header = 10 + (unsigned char)bytes[8] + 256 * (unsigned char)bytes[9];
    return Eigen::Map<const Eigen::Array<T, Eigen::Dynamic, 1>>(reinterpret_cast<const T*>(bytes.data() + header),
                                                                (bytes.size() - header) / sizeof(T));
}

TEST_CASE("fast marching matches scikit-fmm outside the element", "[fmm]")
{
    const std::string fixture = GENERATE("point_source", "ball_variable_speed", "wall_and_dead_zone");
    const int order = GENERATE(1, 2);
    CAPTURE(fixture, order);
    const Eigen::ArrayXd phi = load_npy<double>(fixture + "/phi");
    const pavel::VoxelGrid g{1, Eigen::Vector3i::Zero(), Eigen::Vector3i::Constant(32), phi < 0};
    const Eigen::ArrayXd expected = load_npy<double>(fixture + "/t" + std::to_string(order));
    const Eigen::ArrayXd T = pavel::fast_marching(g, load_npy<double>(fixture + "/speed"),
                                                  load_npy<bool>(fixture + "/mask"), order);

    int reachable = 0, mismatches = 0;
    double error = 0;
    for (Eigen::Index i = 0; i < T.size(); ++i) {
        if (phi[i] < 0) continue;
        mismatches += std::isnan(expected[i]) != std::isinf(T[i]);
        if (std::isnan(expected[i])) continue;
        error = std::max(error, std::abs(T[i] - expected[i]));
        ++reachable;
    }
    CHECK(mismatches == 0);
    CHECK(reachable > 0);
    CHECK(error < 1e-12);
}
