#pragma once

#include "voxel_grid.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <queue>
#include <utility>
#include <vector>

namespace pavel {

namespace detail {

// Port of scikit-fmm's travel-time marcher (v2025.6.23: base_marcher.cpp, distance_marcher.cpp,
// travel_time_marcher.cpp) for a binary phi (element = -1, rest = +1) and dx = 1.
struct Marcher {
    enum Flag : char { Far, Narrow, Frozen, Blocked };
    using Entry = std::pair<double, Eigen::Index>;

    Eigen::Vector3i size;
    const Eigen::Array<bool, Eigen::Dynamic, 1>& inside;
    const Eigen::ArrayXd& speed;
    int order;
    std::vector<Flag> flag;
    Eigen::ArrayXd T;
    // ponytail: lazy deletion (stale entries are skipped) instead of skfmm's decrease-key heap.
    std::priority_queue<Entry, std::vector<Entry>, std::greater<>> heap{};

    Marcher(const VoxelGrid& element, const Eigen::ArrayXd& speed,
            const Eigen::Array<bool, Eigen::Dynamic, 1>& obstacle, int order)
        : size(element.size), inside(element.inside), speed(speed), order(order), flag(size.prod(), Far),
          T(Eigen::ArrayXd::Constant(size.prod(), std::numeric_limits<double>::infinity()))
    {
        for (Eigen::Index i = 0; i < T.size(); ++i)
            if (obstacle[i] || speed[i] < std::numeric_limits<double>::epsilon()) flag[i] = Blocked;
    }

    // Neighbor `step` voxels along `axis`, or -1 if it is off the grid or blocked (skfmm's _getN).
    Eigen::Index open(Eigen::Index i, int axis, int step) const
    {
        const Eigen::Index stride = axis == 0 ? 1 : axis == 1 ? size.x() : Eigen::Index(size.x()) * size.y();
        const Eigen::Index c = i / stride % size[axis] + step;
        if (c < 0 || c >= size[axis] || flag[i + step * stride] == Blocked) return -1;
        return i + step * stride;
    }

    // Eq. 2: Σ_axes max(D⁻T, −D⁺T, 0)² = 1/F², upwinding from the frozen neighbors. Second order
    // uses the longer stencil where two frozen voxels line up and, on a negative discriminant,
    // retries without each axis (bits of `skip`), keeping the smallest time.
    double update(Eigen::Index i, int skip = 0) const
    {
        constexpr double max = std::numeric_limits<double>::max();
        double a = 0, b = 0, c = 0;
        for (int axis = 0; axis < 3; ++axis) {
            if (skip & 1 << axis) continue;
            double t1 = max, t2 = max;
            for (int step : {-1, 1}) {
                const auto n1 = open(i, axis, step);
                if (n1 == -1 || flag[n1] != Frozen || T[n1] >= t1) continue;
                t1 = T[n1];
                // skfmm quirk kept: t2 is not reset when the other side wins t1 without a second voxel.
                // skfmm marches signed times, so a voxel across the interface counts as negative.
                if (const auto n2 = open(i, axis, 2 * step); n2 != -1 && flag[n2] == Frozen && T[n2] <= t1)
                    t2 = inside[n2] != inside[n1] || inside[n2] != inside[i] ? -T[n2] : T[n2];
            }
            if (t1 == max) continue;
            const bool second = order == 2 && t2 < max;
            const double w = second ? 9.0 / 4 : 1, t = second ? (4 * t1 - t2) / 3 : t1;
            a += w;
            b -= 2 * w * t;
            c += w * t * t;
        }
        // larger root of a T² + b T + c = 1/F²
        if (const double det = b * b - 4 * a * (c - 1 / (speed[i] * speed[i])); det >= 0)
            return (-b + std::sqrt(det)) / 2.0 / a;
        if (order == 1) return -b / (2.0 * a); // skfmm: the parabola's minimum as fallback
        double r = std::numeric_limits<double>::infinity();
        for (int axis = 0; axis < 3; ++axis)
            if (!(skip & 1 << axis)) r = std::min(r, update(i, skip | 1 << axis));
        return r;
    }

    void push(Eigen::Index i)
    {
        T[i] = update(i);
        flag[i] = Narrow;
        heap.emplace(T[i], i);
    }

    // Freezes both sides of the interface at the start time 0.5 / (√n F), n = axes crossing it.
    std::vector<Eigen::Index> freeze_boundary()
    {
        std::vector<Eigen::Index> group;
        for (Eigen::Index i = 0; i < T.size(); ++i) {
            if (flag[i] == Blocked) continue;
            int axes = 0;
            for (int axis = 0; axis < 3; ++axis) {
                const auto lo = open(i, axis, -1), hi = open(i, axis, 1);
                axes += (lo != -1 && inside[lo] != inside[i]) || (hi != -1 && inside[hi] != inside[i]);
            }
            if (axes == 0) continue;
            T[i] = std::sqrt(1 / (4.0 * axes)) / speed[i];
            flag[i] = Frozen;
            group.push_back(i);
        }
        return group;
    }

    // As in skfmm, all voxels at exactly the smallest time freeze before any neighbor updates.
    std::vector<Eigen::Index> freeze_next()
    {
        std::vector<Eigen::Index> group;
        const double t = heap.top().first;
        while (!heap.empty() && heap.top().first == t) {
            const auto i = heap.top().second;
            heap.pop();
            if (flag[i] == Narrow && T[i] == t) {
                flag[i] = Frozen;
                group.push_back(i);
            }
        }
        return group;
    }

    // (Re)computes the outside neighbors of newly frozen voxels.
    void expand(const std::vector<Eigen::Index>& group)
    {
        for (const auto i : group)
            for (int axis = 0; axis < 3; ++axis)
                for (int step : {-1, 1}) {
                    const auto n = open(i, axis, step);
                    if (n == -1) continue;
                    if ((flag[n] == Far && !inside[n]) || flag[n] == Narrow) push(n);
                    // second order: the voxel beyond a frozen neighbor now has a longer stencil
                    if (order == 2 && flag[n] == Frozen)
                        if (const auto n2 = open(i, axis, 2 * step); n2 != -1 && flag[n2] == Narrow) push(n2);
                }
    }
};

} // namespace detail

// Arrival times of the front that leaves `element.inside` with per-voxel speed F (§4, Eq. 2),
// in voxel units, matching scikit-fmm's travel_time with phi = -1 inside, +1 outside, dx = 1.
// Obstacles and voxels with F < machine epsilon are never reached. The front starts halfway
// between voxels: a voxel with element neighbors along n axes starts at 0.5 / (√n F).
// Only the outside is marched; T is +inf inside the element and where the front never arrives.
// PAPER: §4/Eq. 2 is first order; scikit-fmm defaults to second order, so `order` switches.
inline Eigen::ArrayXd fast_marching(const VoxelGrid& element, const Eigen::ArrayXd& speed,
                                    const Eigen::Array<bool, Eigen::Dynamic, 1>& obstacle, int order = 1)
{
    detail::Marcher m(element, speed, obstacle, order);
    m.expand(m.freeze_boundary());
    while (!m.heap.empty()) m.expand(m.freeze_next());
    return element.inside.select(std::numeric_limits<double>::infinity(), m.T);
}

} // namespace pavel
