#pragma once

namespace pavel {

// Volume recovery function f_{a,b} (§4, Eq. 1). Maps the distance from a contact surface,
// normalized by the cube root of the element volume, to a front speed.
// The printed Eq. 1 is not C1. This version matches the paper's plot (docs/arxiv/function.png):
// f(0) = 0, f(a) = 1, f((a+b)/2) = 1/2, f(b) = 0, with three quadratics joined C1 at a, (a+b)/2, b.
constexpr double recovery_function(double x, double a, double b)
{
    if (x <= 0 || x > b) return 0;
    if (x <= a) return -x * x / (a * a) + 2 * x / a;
    const double w = (b - a) * (b - a);
    if (x <= (a + b) / 2) return 1 - 2 * (x - a) * (x - a) / w;
    return 2 * (x - b) * (x - b) / w;
}

} // namespace pavel
