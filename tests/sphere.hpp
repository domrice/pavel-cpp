#pragma once

#include <igl/upsample.h>

// Unit sphere: an octahedron upsampled 4 times (2048 faces), projected onto the sphere.
inline void unit_sphere(Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    V.resize(6, 3);
    V << 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1;
    F.resize(8, 3);
    F << 0, 2, 4, 2, 1, 4, 1, 3, 4, 3, 0, 4, 2, 0, 5, 1, 2, 5, 3, 1, 5, 0, 3, 5;
    igl::upsample(Eigen::MatrixXd(V), Eigen::MatrixXi(F), V, F, 4);
    V.rowwise().normalize();
}
