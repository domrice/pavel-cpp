#pragma once

#include <Eigen/Core>

namespace pavel {

using Mask = Eigen::Array<bool, Eigen::Dynamic, 1>;

// A dense block of the global lattice (§4): every grid uses the same voxel size h and the
// world axes, and the voxel with lattice index n has its center at n * h. A grid is then just an
// integer offset plus extents, so Booleans between grids are index arithmetic.
struct VoxelGrid {
    double h;
    Eigen::Vector3i offset; // lattice index of the grid's voxel (0, 0, 0)
    Eigen::Vector3i size;
    Mask inside; // x fastest, as igl::marching_cubes expects

    // Lattice index of the voxel stored at i.
    Eigen::Vector3i coord(Eigen::Index i) const
    {
        return offset + Eigen::Vector3i(i % size.x(), i / size.x() % size.y(), i / (size.x() * size.y()));
    }

    // Storage index of lattice voxel p, or -1 if p is off the grid.
    Eigen::Index index(const Eigen::Vector3i& p) const
    {
        const Eigen::Vector3i q = p - offset;
        if ((q.array() < 0).any() || (q.array() >= size.array()).any()) return -1;
        return q.x() + Eigen::Index(size.x()) * (q.y() + Eigen::Index(size.y()) * q.z());
    }

    // Whether lattice voxel p is inside (false off the grid).
    bool at(const Eigen::Vector3i& p) const
    {
        const Eigen::Index i = index(p);
        return i != -1 && inside[i];
    }

    // other.inside at this grid's voxels shifted by s (false where other does not reach).
    Mask sample(const VoxelGrid& other, const Eigen::Vector3i& s = Eigen::Vector3i::Zero()) const
    {
        Mask m(size.prod());
        for (Eigen::Index i = 0; i < m.size(); ++i) m[i] = other.at(coord(i) + s);
        return m;
    }

    // World positions of all voxel centers, in storage order.
    Eigen::MatrixXd centers() const
    {
        Eigen::MatrixXd C(size.prod(), 3);
        for (int i = 0; i < C.rows(); ++i) C.row(i) = h * coord(i).cast<double>().transpose();
        return C;
    }
};

} // namespace pavel
