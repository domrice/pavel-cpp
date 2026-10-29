#pragma once

#include "fast_marching.hpp"
#include "recovery_function.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

namespace pavel {

// §4 defaults: a = 0.1, b = 0.3 for every paper result except Fig. 5.
struct DeformParams {
    double a = 0.1, b = 0.3;
    int order = 1; // of the fast marching
};

// An element whose grid intersects element i's, so the two compete for voxels.
struct Neighbor {
    int j;
    int lost = 0; // voxels i lost that j also covered; i and j are in contact if > 0
};

// The elements after overlap resolution (§4), D′.
struct Overlaps {
    std::vector<VoxelGrid> resolved;
    std::vector<std::vector<Neighbor>> neighbors;
};

// Removes element voxels inside the base, and voxels shared by several elements from all but the one with the
// nearest center (ties go to the lower index).
// PAPER: §4, the element center is the centroid of its voxels.
inline Overlaps resolve_overlaps(const std::vector<VoxelGrid>& elements, const VoxelGrid& base)
{
    const int n = elements.size();
    Overlaps o{elements, std::vector<std::vector<Neighbor>>(n)};
    // ponytail: O(n²) grid box tests; a spatial hash once n is in the thousands.
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j) {
            const VoxelGrid &gi = elements[i], &gj = elements[j];
            if ((gi.offset.array() < (gj.offset + gj.size).array()).all() &&
                (gj.offset.array() < (gi.offset + gi.size).array()).all()) {
                o.neighbors[i].push_back({j});
                o.neighbors[j].push_back({i});
            }
        }
    std::vector<Eigen::RowVector3d> center(n);
    for (int i = 0; i < n; ++i)
        center[i] = elements[i].inside.cast<double>().matrix().transpose() * elements[i].centers() /
                    elements[i].inside.count();

    for (int i = 0; i < n; ++i) {
        const VoxelGrid& g = elements[i];
        for (Eigen::Index v = 0; v < g.inside.size(); ++v) {
            if (!g.inside[v]) continue;
            const Eigen::Vector3i p = g.coord(v);
            const Eigen::RowVector3d x = g.h * p.cast<double>().transpose();
            const std::pair mine{(x - center[i]).squaredNorm(), i};
            const bool removed = base.at(p) || std::ranges::any_of(o.neighbors[i], [&](const Neighbor& nb) {
                return elements[nb.j].at(p) && std::pair{(x - center[nb.j]).squaredNorm(), nb.j} < mine;
            });
            if (!removed) continue;
            o.resolved[i].inside[v] = false;
            for (auto& nb : o.neighbors[i]) nb.lost += elements[nb.j].at(p);
        }
    }
    return o;
}

// Whether lattice voxel p belongs to a neighbor of element i.
inline bool taken(const Overlaps& o, int i, const Eigen::Vector3i& p)
{
    return std::ranges::any_of(o.neighbors[i], [&](const Neighbor& nb) { return o.resolved[nb.j].at(p); });
}

// Phase A for element i: the speed F = Σ_c ω_c f_{a,b}(r_c / d) (Eq. 3) and the arrival times of the front leaving
// the element (Eq. 2), in voxel units. The base and the other elements are obstacles.
// PAPER: §4, ω_c is the voxel count i lost in its overlap with c (its scale does not matter in round-robin phase B),
// r_c is the distance to i's voxels touching c, and the base is not a contact.
inline Eigen::ArrayXd arrival_times(const Overlaps& o, int i, const VoxelGrid& base, double d, const DeformParams& p)
{
    const VoxelGrid& g = o.resolved[i];
    const Eigen::Index n = g.inside.size();
    Eigen::ArrayXd F = Eigen::ArrayXd::Zero(n);
    Mask obstacle = g.sample(base);
    for (auto [j, lost] : o.neighbors[i]) {
        const VoxelGrid& c = o.resolved[j];
        obstacle = obstacle || g.sample(c);
        if (lost == 0) continue;
        Mask touch = Mask::Constant(n, false);
        for (int axis = 0; axis < 3; ++axis)
            for (int step : {-1, 1}) touch = touch || g.sample(c, step * Eigen::Vector3i::Unit(axis));
        const VoxelGrid contact{g.h, g.offset, g.size, g.inside && touch};
        const Eigen::ArrayXd r = fast_marching(contact, Eigen::ArrayXd::Ones(n), Mask::Constant(n, false), p.order);
        for (Eigen::Index v = 0; v < n; ++v) F[v] += lost * recovery_function(r[v] / d, p.a, p.b);
    }
    return fast_marching(g, F, obstacle, p.order);
}

// Queues the reachable 6-neighbors of voxel v that are not queued yet. A queued voxel's T becomes +inf, so T
// doubles as the queued mark.
inline void push_neighbors(const VoxelGrid& g, Eigen::ArrayXd& T, MinHeap& front, Eigen::Index v)
{
    for (int axis = 0; axis < 3; ++axis)
        for (int step : {-1, 1}) {
            const Eigen::Index u = g.index(g.coord(v) + step * Eigen::Vector3i::Unit(axis));
            if (u == -1 || std::isinf(T[u])) continue;
            front.emplace(std::exchange(T[u], std::numeric_limits<double>::infinity()), u);
        }
}

// Volume recovery (§4): resolves the overlaps, then grows each element back to its original voxel count into free
// voxels, in order of its arrival times. Phase A computes all arrival times concurrently, phase B assigns voxels
// sequentially. An element whose front runs out of reachable voxels stays short (§4: "theoretically possible").
// PAPER: §4, phase B is round-robin (arXiv text): each element in turn claims one voxel.
inline std::vector<VoxelGrid> deform(const std::vector<VoxelGrid>& elements, const VoxelGrid& base,
                                     const DeformParams& p = {})
{
    Overlaps o = resolve_overlaps(elements, base);
    const int n = elements.size();

    std::vector<Eigen::ArrayXd> T(n);
    std::atomic<int> next = 0;
    std::vector<std::thread> workers(std::max(1u, std::thread::hardware_concurrency()));
    for (auto& w : workers)
        w = std::thread([&] {
            for (int i; (i = next++) < n;) T[i] = arrival_times(o, i, base, std::cbrt(elements[i].inside.count()), p);
        });
    for (auto& w : workers) w.join();

    std::vector<MinHeap> fronts(n);
    std::vector<Eigen::Index> deficit(n);
    for (int i = 0; i < n; ++i) {
        const VoxelGrid& g = o.resolved[i];
        deficit[i] = elements[i].inside.count() - g.inside.count();
        for (Eigen::Index v = 0; v < g.inside.size(); ++v)
            if (g.inside[v]) push_neighbors(g, T[i], fronts[i], v);
    }
    for (bool claimed = true; claimed;) {
        claimed = false;
        for (int i = 0; i < n; ++i) {
            VoxelGrid& g = o.resolved[i];
            MinHeap& front = fronts[i];
            if (deficit[i] == 0) continue;
            while (!front.empty() && taken(o, i, g.coord(front.top().second))) front.pop();
            if (front.empty()) continue;
            const Eigen::Index v = front.top().second;
            front.pop();
            g.inside[v] = true;
            --deficit[i];
            claimed = true;
            push_neighbors(g, T[i], front, v);
        }
    }
    return std::move(o.resolved);
}

} // namespace pavel
