#pragma once

#include <Eigen/Core>
#include <geogram/basic/command_line_args.h>
#include <geogram/basic/common.h>
#include <geogram/mesh/mesh.h>

namespace pavel {

// Triangle mesh V, F into M. Also initializes Geogram (idempotent), since every use of it starts here.
inline void to_geogram(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, GEO::Mesh& M)
{
    GEO::initialize();
    GEO::CmdLine::import_arg_group("standard");
    GEO::CmdLine::import_arg_group("algo");
    const Eigen::Matrix<double, 3, Eigen::Dynamic> Vt = V.transpose();
    M.vertices.assign_points(Vt.data(), 3, GEO::index_t(V.rows()));
    GEO::vector<GEO::index_t> T(3 * F.rows());
    Eigen::Map<Eigen::Matrix<GEO::index_t, 3, Eigen::Dynamic>>(T.data(), 3, F.rows()) =
        F.transpose().cast<GEO::index_t>();
    M.facets.assign_triangle_mesh(T, true);
    M.facets.connect();
}

} // namespace pavel
