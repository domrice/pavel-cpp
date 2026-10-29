#include "seeding.hpp"
#include "sphere.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <sstream>

using Catch::Matchers::WithinRel;

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

// Smallest distance between two seeds.
static double min_spacing(const std::vector<pavel::Seed>& seeds)
{
    double d = INFINITY;
    for (size_t i = 0; i < seeds.size(); ++i)
        for (size_t j = 0; j < i; ++j) d = std::min(d, (seeds[i].position - seeds[j].position).norm());
    return d;
}

TEST_CASE("isotropic seeds cover a sphere evenly", "[seeding]")
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    unit_sphere(V, F);
    const double fp = pavel::footprint_area(V, F, 0.02); // the equator disk
    CHECK_THAT(fp, WithinRel(std::numbers::pi, 0.02));

    // Unit-sphere elements on a radius-10 base at coverage 1.4: 1.4 · 400π / π = 560 seeds.
    const Eigen::MatrixXd base = 10 * V;
    const auto seeds = pavel::isotropic_seeds(base, F, fp, {});
    CHECK_THAT(double(seeds.size()), WithinRel(560, 0.03));
    for (const auto& s : seeds) {
        CHECK_THAT(s.position.norm(), WithinRel(10, 0.01));
        CHECK(s.normal.dot(s.position.normalized()) > 0.999);
    }

    // Converged, the seeds stay near the spacing of a hexagonal packing; with no iterations they are random and clump.
    const double hex = std::sqrt(2 * 400 * std::numbers::pi / (std::sqrt(3) * seeds.size()));
    const double cvt = min_spacing(seeds), random = min_spacing(pavel::isotropic_seeds(base, F, fp, {.lloyd = 0, .newton = 0}));
    INFO("hexagonal " << hex << " cvt " << cvt << " random " << random);
    CHECK(cvt > 0.75 * hex);
    CHECK(random < 0.5 * cvt);
}
