#pragma once

#include <Eigen/Core>

namespace pavel {

// A dense block of the global lattice (§4): every grid uses the same voxel size h and the
// world axes, and the voxel with lattice index n has its center at n * h. A grid is then just an
// integer offset plus extents, so Booleans between grids are index arithmetic.
struct VoxelGrid {
    double h;
    Eigen::Vector3i offset; // lattice index of the grid's voxel (0, 0, 0)
    Eigen::Vector3i size;
    Eigen::Array<bool, Eigen::Dynamic, 1> inside; // x fastest, as igl::marching_cubes expects

    // World positions of all voxel centers, in storage order.
    Eigen::MatrixXd centers() const
    {
        Eigen::MatrixXd C(size.prod(), 3);
        for (int i = 0; i < C.rows(); ++i)
            C.row(i) << i % size.x(), i / size.x() % size.y(), i / (size.x() * size.y());
        return h * (C.rowwise() + offset.cast<double>().transpose());
    }
};

} // namespace pavel
