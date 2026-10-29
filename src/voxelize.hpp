#pragma once

#include "voxel_grid.hpp"

#include <Eigen/Geometry>
#include <igl/fast_winding_number.h>
#include <igl/marching_cubes.h>

namespace pavel {

// PAPER: §4, the grid is 50% larger than the element: its bounding box scaled by 1.5 about
// the center.
inline Eigen::AlignedBox3d element_box(const Eigen::MatrixXd& V)
{
    const Eigen::AlignedBox3d b(V.colwise().minCoeff(), V.colwise().maxCoeff());
    return {b.center() - 0.75 * b.sizes(), b.center() + 0.75 * b.sizes()};
}

// Voxelizes a closed mesh on the lattice of voxel size h, over the voxels needed to cover box
// (rounded outward to whole voxels). A voxel is inside when the fast winding number at its
// center exceeds 1/2.
inline VoxelGrid voxelize(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double h, const Eigen::AlignedBox3d& box)
{
    const Eigen::Vector3i lo = (box.min() / h).array().floor().cast<int>().matrix();
    const Eigen::Vector3i hi = (box.max() / h).array().ceil().cast<int>().matrix();
    VoxelGrid g{h, lo, hi - lo + Eigen::Vector3i::Ones(), {}};
    Eigen::VectorXd W;
    igl::fast_winding_number(V, F, g.centers(), W);
    g.inside = W.array() > 0.5;
    return g;
}

// Surface of the inside voxels: marching cubes on the 0/1 field at isovalue 1/2, so the surface
// passes halfway between inside and outside voxel centers. The field is negated (inside = -1)
// because igl orients faces toward larger values, and we want outward normals.
inline void marching_cubes(const VoxelGrid& g, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    igl::marching_cubes(-g.inside.cast<double>().matrix(), g.centers(), g.size.x(), g.size.y(), g.size.z(),
                        -0.5, V, F);
}

} // namespace pavel
