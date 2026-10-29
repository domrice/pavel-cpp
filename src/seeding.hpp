#pragma once

#include "mesh_convert.hpp"
#include "voxelize.hpp"

#include <Eigen/Geometry>
#include <geogram/voronoi/CVT.h>
#include <igl/barycentric_coordinates.h>
#include <igl/doublearea.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>

#include <cmath>
#include <istream>
#include <numbers>
#include <stdexcept>
#include <vector>

namespace pavel {

// Where one element goes (§5): its frame origin at position, its +z along normal, turned by angle (degrees) about
// the normal, and scaled uniformly. element indexes the list of element meshes.
struct Seed {
    Eigen::Vector3d position, normal;
    double angle;
    int element;
    double scale;
};

// Manual seeds, whitespace separated, one per line: x y z nx ny nz angle element scale.
inline std::vector<Seed> read_seeds(std::istream& in)
{
    std::vector<Seed> seeds;
    while (!(in >> std::ws).eof()) {
        Seed s;
        if (!(in >> s.position.x() >> s.position.y() >> s.position.z() >> s.normal.x() >> s.normal.y() >>
              s.normal.z() >> s.angle >> s.element >> s.scale))
            throw std::runtime_error("malformed seeds file");
        seeds.push_back(s);
    }
    return seeds;
}

// Element vertices V, given in the element frame (+z the symmetry axis, z = 0 the footprint plane), placed at s.
// Angle 0 is the minimal rotation taking +z to the normal.
inline Eigen::MatrixXd place(const Eigen::MatrixXd& V, const Seed& s)
{
    const Eigen::Matrix3d R = (Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), s.normal) *
                               Eigen::Quaterniond(Eigen::AngleAxisd(s.angle * std::numbers::pi / 180,
                                                                    Eigen::Vector3d::UnitZ())))
                                  .toRotationMatrix();
    return (s.scale * V * R.transpose()).rowwise() + s.position.transpose();
}

// PAPER: §4, the footprint of an element (V, F in its frame) is its z-section of maximal area, measured on its voxels
// of size h.
inline double footprint_area(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double h)
{
    const auto g = voxelize(V, F, h, element_box(V));
    const auto layers = g.inside.cast<int>().reshaped(Eigen::Index(g.size.x()) * g.size.y(), g.size.z());
    return h * h * layers.colwise().sum().maxCoeff();
}

// PAPER: §4, coverage about 1.4 (1.2 barely deforms, 1.6 distorts). The iteration counts are GEO::remesh_smooth's.
struct IsotropicParams {
    double coverage = 1.4;
    int lloyd = 5, newton = 30;
};

// Isotropic seeding (§5) on the base V, F: coverage = N · footprint / base area fixes the number N of seeds, which
// Geogram's CVT (Yan et al. 2009) spreads over the surface with lloyd then newton iterations; fewer give a less
// regular, "hand-made" look. Each seed is projected onto the base, its normal interpolates the base's vertex normals,
// and it places element 0 unturned and unscaled.
inline std::vector<Seed> isotropic_seeds(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double footprint,
                                         const IsotropicParams& p)
{
    Eigen::VectorXd A;
    igl::doublearea(V, F, A);
    const auto n = GEO::index_t(std::lround(p.coverage * A.sum() / 2 / footprint));
    GEO::Mesh M;
    to_geogram(V, F, M);
    GEO::CentroidalVoronoiTesselation cvt(&M);
    cvt.compute_initial_sampling(n);
    cvt.Lloyd_iterations(p.lloyd);
    if (p.newton > 0) cvt.Newton_iterations(p.newton); // as in GEO::remesh_smooth

    const Eigen::MatrixXd P = Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(
        cvt.embedding(0), n, 3);
    Eigen::MatrixXd C, B, N;
    Eigen::VectorXd D;
    Eigen::VectorXi I;
    igl::point_mesh_squared_distance(P, V, F, D, I, C);
    using Eigen::placeholders::all;
    igl::barycentric_coordinates(C, V(F(I, 0), all), V(F(I, 1), all), V(F(I, 2), all), B);
    igl::per_vertex_normals(V, F, N);
    std::vector<Seed> seeds;
    for (GEO::index_t i = 0; i < n; ++i) {
        const Eigen::Vector3d normal = B.row(i) * N(F.row(I[i]), all);
        seeds.push_back({C.row(i), normal.normalized(), 0, 0, 1});
    }
    return seeds;
}

} // namespace pavel
