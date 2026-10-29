#include "voxelize.hpp"

#include <igl/centroid.h>
#include <igl/upsample.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <numbers>

using Catch::Matchers::WithinRel;

TEST_CASE("unit sphere survives voxelize -> marching cubes", "[voxelize]")
{
    Eigen::MatrixXd V(6, 3);
    V << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    Eigen::MatrixXi F(8, 3);
    F << 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5;
    igl::upsample(Eigen::MatrixXd(V), Eigen::MatrixXi(F), V, F, 4);
    V.rowwise().normalize();

    const double h = 0.02; // diagonal of the 1.5x grid: about 225 voxels, like the paper's
    const auto g = pavel::voxelize(V, F, h, pavel::element_box(V));
    Eigen::MatrixXd MV;
    Eigen::MatrixXi MF;
    pavel::marching_cubes(g, MV, MF);

    Eigen::Vector3d centroid;
    double mesh, mc;
    igl::centroid(V, F, centroid, mesh);
    igl::centroid(MV, MF, centroid, mc);
    const double sphere = 4 * std::numbers::pi / 3, voxels = g.inside.count() * h * h * h;
    INFO("sphere " << sphere << " mesh " << mesh << " voxels " << voxels << " mc " << mc);
    CHECK((g.size.cast<double>().array() * h > 3).all()); // 1.5x the diameter
    CHECK_THAT(voxels, WithinRel(mesh, 0.01));
    CHECK_THAT(mc, WithinRel(mesh, 0.01));
    CHECK((MV.rowwise().norm().array() - 1).abs().maxCoeff() < h);
}
