#include "recovery_function.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using Catch::Matchers::WithinAbs;
using pavel::recovery_function;

constexpr double a = 0.1, b = 0.3;

TEST_CASE("recovery function hits the control points of the plot", "[recovery]")
{
    CHECK(recovery_function(0, a, b) == 0);
    CHECK_THAT(recovery_function(a, a, b), WithinAbs(1, 1e-12));
    CHECK_THAT(recovery_function((a + b) / 2, a, b), WithinAbs(0.5, 1e-12));
    CHECK(recovery_function(b, a, b) == 0);
    CHECK(recovery_function(1, a, b) == 0);
}

TEST_CASE("recovery function is C1 at the joins", "[recovery]")
{
    constexpr double h = 1e-8;
    for (double x : {a, (a + b) / 2, b}) {
        const double left = (recovery_function(x, a, b) - recovery_function(x - h, a, b)) / h;
        const double right = (recovery_function(x + h, a, b) - recovery_function(x, a, b)) / h;
        CHECK_THAT(left, WithinAbs(right, 1e-5));
    }
}
