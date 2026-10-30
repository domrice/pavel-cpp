#include "output.hpp"
#include "sphere.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinRel;

TEST_CASE("remeshing two overlapping spheres keeps them and the point count", "[output]")
{
    Eigen::MatrixXd S, RV;
    Eigen::MatrixXi SF, RF;
    unit_sphere(S, SF);
    const Eigen::RowVector3d shift(1, 0, 0);
    const std::vector<Eigen::MatrixXd> V{S, S.rowwise() + shift};
    const std::vector<Eigen::MatrixXi> F{SF, SF};

    // A seed whose restricted Voronoi cell is split (here, across the crossing) gives one vertex per piece.
    pavel::remesh(V, F, 1000, RV, RF);
    CHECK_THAT(double(RV.rows()), WithinRel(1000, 0.1));
    REQUIRE(RF.rows() > 0);
    CHECK(RF.maxCoeff() < RV.rows());
    for (Eigen::Index i = 0; i < RV.rows(); ++i) {
        const double d = std::min(std::abs(RV.row(i).norm() - 1), std::abs((RV.row(i) - shift).norm() - 1));
        CHECK(d < 0.01);
    }

    pavel::remesh(V, F, 0, RV, RF);
    CHECK_THAT(double(RV.rows()), WithinRel(2.0 * S.rows(), 0.05));
}

TEST_CASE("patches split the seeds by the base quadrant", "[output]")
{
    // Only the base's bounding box counts: centered at (1, 1, 0), so the cut planes are x = 1 and y = 1. The last seed
    // lies on both.
    const Eigen::MatrixXd BV = (Eigen::MatrixXd(2, 3) << -9, -9, -10, 11, 11, 10).finished();
    std::vector<pavel::Seed> seeds;
    for (const auto& p : {Eigen::Vector3d(-2, -2, 0), Eigen::Vector3d(4, -2, 0), Eigen::Vector3d(-2, 4, 0),
                          Eigen::Vector3d(4, 4, 0), Eigen::Vector3d(1, 1, 0)})
        seeds.push_back({p, Eigen::Vector3d::UnitZ(), 0, 0, 1});
    CHECK(pavel::patches(BV, seeds) == std::vector<int>{0, 1, 2, 3, 3});
}
