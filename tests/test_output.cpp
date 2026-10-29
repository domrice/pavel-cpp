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
