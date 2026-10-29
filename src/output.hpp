#pragma once

#include "mesh_convert.hpp"

#include <geogram/mesh/mesh_remesh.h>
#include <igl/combine.h>

namespace pavel {

// PAPER: §6, the elements' marching-cubes meshes and the base, V[k], F[k], are merged into one model (a plain
// concatenation: faces shared by touching elements and element bottoms inside the base stay), and CVT-remeshed
// (Yan et al. 2009) with n points, which is also the smoothing step. n = 0 keeps the merged vertex count.
inline void remesh(const std::vector<Eigen::MatrixXd>& V, const std::vector<Eigen::MatrixXi>& F, unsigned n,
                   Eigen::MatrixXd& RV, Eigen::MatrixXi& RF)
{
    Eigen::MatrixXd CV;
    Eigen::MatrixXi CF;
    igl::combine(V, F, CV, CF);
    GEO::Mesh in, out;
    to_geogram(CV, CF, in);
    GEO::remesh_smooth(in, out, n);
    from_geogram(out, RV, RF);
}

} // namespace pavel
