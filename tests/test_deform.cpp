#include "deform.hpp"
#include "voxelize.hpp"

#include <igl/writeOBJ.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>

#include <filesystem>
#include <format>

using pavel::VoxelGrid;

constexpr double h = 0.05;

// Ball of radius r around c in a grid 50% larger than the ball.
VoxelGrid ball(const Eigen::Vector3d& c, double r)
{
    VoxelGrid g{h, ((c.array() - 1.5 * r) / h).floor().cast<int>(), Eigen::Vector3i::Constant(3 * r / h + 2), {}};
    g.inside = (g.centers().rowwise() - c.transpose()).rowwise().norm().array() < r;
    return g;
}

// Two unit balls overlapping each other and sunk 0.1 into the base half-space z <= 0, as in Fig. 5.
const std::vector<VoxelGrid> elements{ball({-0.8, 0, 0.9}, 1), ball({0.8, 0, 0.9}, 1)};
const VoxelGrid base{h, {-50, -50, -20}, {100, 100, 21}, pavel::Mask::Constant(100 * 100 * 21, true)};

TEST_CASE("overlap resolution leaves no voxel shared", "[deform]")
{
    const auto o = pavel::resolve_overlaps(elements, base);
    for (int i : {0, 1}) {
        REQUIRE(o.neighbors[i].size() == 1);
        CHECK(o.neighbors[i][0].j == 1 - i);
        CHECK(o.neighbors[i][0].lost > 0);
        CHECK(o.resolved[i].inside.count() < elements[i].inside.count() - o.neighbors[i][0].lost); // base took some too
    }
}

TEST_CASE("deformation recovers the exact volume without sharing voxels", "[deform]")
{
    const auto deformed = pavel::deform(elements, base);
    for (int i : {0, 1}) {
        CHECK(deformed[i].inside.count() == elements[i].inside.count());
        int shared = 0;
        for (Eigen::Index v = 0; v < deformed[i].inside.size(); ++v)
            if (deformed[i].inside[v])
                shared += base.at(deformed[i].coord(v)) || deformed[1 - i].at(deformed[i].coord(v));
        CHECK(shared == 0);
    }
}

TEST_CASE("Fig. 5 meshes for a visual check", "[.][visual]")
{
    const auto [a, b] = GENERATE(std::pair{0.05, 0.15}, std::pair{0.1, 0.3}, std::pair{0.2, 0.6});
    CAPTURE(a, b);
    const auto deformed = pavel::deform(elements, base, {a, b});
    std::filesystem::create_directories(PAVEL_OUT);
    for (int i : {0, 1}) {
        CHECK(deformed[i].inside.count() == elements[i].inside.count());
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        pavel::marching_cubes(deformed[i], V, F);
        igl::writeOBJ(std::format("{}/two_balls_a{}_b{}_{}.obj", PAVEL_OUT, a, b, i), V, F);
    }
}
