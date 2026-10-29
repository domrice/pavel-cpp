#include "seeding.hpp"

#include <catch2/catch_test_macros.hpp>

#include <sstream>

TEST_CASE("seeds are read and place elements", "[seeding]")
{
    std::istringstream in("1 2 3  1 0 0  90 1 2\n0 0 0  0 0 1  0 0 1\n");
    const auto seeds = pavel::read_seeds(in);
    REQUIRE(seeds.size() == 2);
    CHECK(seeds[0].element == 1);

    // +z goes along the normal +x; the 90° turn takes +x to +y, which the z → x alignment (about y) keeps.
    const Eigen::MatrixXd V = pavel::place(Eigen::MatrixXd{{0, 0, 1}, {1, 0, 0}}, seeds[0]);
    CHECK((V - Eigen::MatrixXd{{3, 2, 3}, {1, 4, 3}}).norm() < 1e-12);

    std::istringstream bad("1 2 3 x"), truncated("1 2 3\n");
    CHECK_THROWS(pavel::read_seeds(bad));
    CHECK_THROWS(pavel::read_seeds(truncated));
}
