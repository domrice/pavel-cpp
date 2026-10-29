#include "seeding.hpp"
#include "sphere.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <igl/triangulated_grid.h>

#include <sstream>

using Catch::Matchers::WithinRel;

TEST_CASE("seeds are read and place elements", "[seeding]")
{
    std::istringstream in("1 2 3  1 0 0  90 1 2\n0 0 0  0 0 1  0 0 1\n");
    const auto seeds = pavel::read_seeds(in);
    REQUIRE(seeds.size() == 2);
    CHECK(seeds[0].element == 1);
    std::ostringstream out;
    pavel::write_seeds(out, seeds);
    CHECK(out.str() == "1 2 3 1 0 0 90 1 2\n0 0 0 0 0 1 0 0 1\n");

    // +z goes along the normal +x; the 90° turn takes +x to +y, which the z → x alignment (about y) keeps.
    const Eigen::MatrixXd V = pavel::place(Eigen::MatrixXd{{0, 0, 1}, {1, 0, 0}}, seeds[0]);
    CHECK((V - Eigen::MatrixXd{{3, 2, 3}, {1, 4, 3}}).norm() < 1e-12);

    // angle_toward inverts place() on +x, also off the axes.
    const Eigen::Vector3d n = Eigen::Vector3d(1, -2, 0.5).normalized();
    const Eigen::Vector3d t = pavel::place(Eigen::MatrixXd{{1, 0, 0}}, {Eigen::Vector3d::Zero(), n, 130, 0, 1}).row(0);
    CHECK_THAT(pavel::angle_toward(n, t), WithinRel(130, 1e-12));

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
    const auto [fp, height] = pavel::footprint({V}, {F}, 0.02); // the equator disk
    CHECK_THAT(fp, WithinRel(std::numbers::pi, 0.02));
    CHECK(std::abs(height) <= 0.02);

    // With a radius-2 sphere centered at z = 1: the mean of π and 4π, and of heights 0 and 1.
    const Eigen::MatrixXd W = (2 * V).rowwise() + Eigen::RowVector3d(0, 0, 1);
    const auto mean = pavel::footprint({V, W}, {F, F}, 0.05);
    CHECK_THAT(mean.area, WithinRel(2.5 * std::numbers::pi, 0.05));
    CHECK(std::abs(mean.height - 0.5) <= 0.05);

    // Unit-sphere elements on a radius-10 base at coverage 1.4: 1.4 · 400π / π = 560 seeds.
    const Eigen::MatrixXd base = 10 * V;
    const auto seeds = pavel::isotropic_seeds(base, F, base, F, fp, {});
    CHECK_THAT(double(seeds.size()), WithinRel(560, 0.03));
    for (const auto& s : seeds) {
        CHECK_THAT(s.position.norm(), WithinRel(10, 0.01));
        CHECK(s.normal.dot(s.position.normalized()) > 0.999);
    }

    // Converged, the seeds stay near the spacing of a hexagonal packing; with no iterations they are random and clump.
    const double hex = std::sqrt(2 * 400 * std::numbers::pi / (std::sqrt(3) * seeds.size()));
    const double cvt = min_spacing(seeds),
                 random = min_spacing(pavel::isotropic_seeds(base, F, base, F, fp, {.lloyd = 0, .newton = 0}));
    INFO("hexagonal " << hex << " cvt " << cvt << " random " << random);
    CHECK(cvt > 0.75 * hex);
    CHECK(random < 0.5 * cvt);
}

TEST_CASE("the offset surface of a sphere is a larger sphere", "[seeding]")
{
    Eigen::MatrixXd V, OV;
    Eigen::MatrixXi F, OF;
    unit_sphere(V, F);
    pavel::offset_surface(V, F, 0.5, 0.05, OV, OF);
    REQUIRE(OV.rows() > 0);
    CHECK((OV.rowwise().norm().array() - 1.5).abs().maxCoeff() < 0.01);
}

TEST_CASE("stripe seeds on a square form a grid alternating along the main stripes", "[seeding]")
{
    Eigen::MatrixXd G;
    Eigen::MatrixXi F;
    igl::triangulated_grid(41, 41, G, F);
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(G.rows(), 3);
    V.leftCols(2) = 10 * G;
    std::vector<int> parity;
    const auto seeds = pavel::stripe_seeds(V, F, 1, 1, 90, parity);
    CHECK_THAT(double(seeds.size()), WithinRel(100, 0.25));
    CHECK(min_spacing(seeds) > 0.8);

    // Neighbours one spacing apart differ in parity along the main stripe (+x of the element) and match across it.
    int along = 0, across = 0, wrong = 0;
    for (size_t i = 0; i < seeds.size(); ++i) {
        const double a = seeds[i].angle * std::numbers::pi / 180;
        const Eigen::Vector3d t(std::cos(a), std::sin(a), 0);
        for (size_t j = 0; j < seeds.size(); ++j) {
            const Eigen::Vector3d e = seeds[j].position - seeds[i].position;
            if (i == j || e.norm() > 1.2) continue;
            const double c = std::abs(e.normalized().dot(t));
            if (c > 0.9) ++along, wrong += parity[i] == parity[j];
            if (c < 0.1) ++across, wrong += parity[i] != parity[j];
        }
    }
    CHECK(along > 100);
    CHECK(across > 100);
    CHECK(wrong == 0);
}

TEST_CASE("perturbation alternates and jitters", "[seeding]")
{
    std::vector<pavel::Seed> seeds(4, {Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitZ(), 0, 0, 1});
    pavel::perturb(seeds, {0, 1, 0, 1}, 2,
                   {.rotation = pavel::Variation::alternate, .element = pavel::Variation::alternate, .scale_jitter = 0.1});
    for (int i = 0; i < 4; ++i) {
        CHECK(seeds[i].angle == 180 * (i % 2));
        CHECK(seeds[i].element == i % 2);
        CHECK(std::abs(seeds[i].scale - 1) <= 0.1);
    }
}
