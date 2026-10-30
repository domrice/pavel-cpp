#pragma once

#include "deform.hpp"
#include "mesh_convert.hpp"
#include "seeding.hpp"

#include <geogram/mesh/mesh_remesh.h>
#include <igl/combine.h>

#include <array>
#include <iostream>

namespace pavel {

// The pipeline from seeds to meshes: the elements EV[s.element], EF[s.element] placed at each seed and voxelized at
// voxel size h, with the base BV, BF voxelized over each element's grid, then deformed (§4), and each one's
// marching-cubes mesh in V[i], F[i]. Elements that end short are reported on stdout.
// ponytail: the base's winding-number tree is rebuilt per element; precompute it once if that shows up in a profile.
inline void deformed_meshes(const Eigen::MatrixXd& BV, const Eigen::MatrixXi& BF, const std::vector<Eigen::MatrixXd>& EV,
                            const std::vector<Eigen::MatrixXi>& EF, const std::vector<Seed>& seeds, double h,
                            const DeformParams& p, std::vector<Eigen::MatrixXd>& V, std::vector<Eigen::MatrixXi>& F)
{
    std::vector<VoxelGrid> elements;
    std::vector<Mask> base;
    for (const auto& s : seeds) {
        const Eigen::MatrixXd P = place(EV.at(s.element), s);
        const auto box = element_box(P);
        elements.push_back(voxelize(P, EF.at(s.element), h, box));
        base.push_back(voxelize(BV, BF, h, box).inside);
    }
    const auto deformed = deform(elements, base, p);
    V.resize(deformed.size());
    F.resize(deformed.size());
    for (size_t i = 0; i < deformed.size(); ++i) {
        if (deformed[i].inside.count() < elements[i].inside.count())
            std::cout << "element " << i << " recovered " << deformed[i].inside.count() << " of "
                      << elements[i].inside.count() << " voxels\n";
        marching_cubes(deformed[i], V[i], F[i]);
    }
}

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

// PAPER: §6, printable shells: the base BV is cut by two perpendicular planes through its symmetry axis, taken as the
// z axis through its bounding-box center c, into the patches k = [x ≥ c.x] + 2 [y ≥ c.y]. Returns each seed's patch.
inline std::vector<int> patches(const Eigen::MatrixXd& BV, const std::vector<Seed>& seeds)
{
    const Eigen::Vector3d c = (BV.colwise().minCoeff() + BV.colwise().maxCoeff()).transpose() / 2;
    std::vector<int> k;
    for (const auto& s : seeds) k.push_back((s.position.x() >= c.x()) + 2 * (s.position.y() >= c.y()));
    return k;
}

// Each element (mesh V[i], F[i]) joins the patch holding its seed, and each patch's elements, without the base, are
// remeshed (keeping their vertex count) into one shell SV[k], SF[k]. A patch without elements leaves its shell empty.
inline void shells(const Eigen::MatrixXd& BV, const std::vector<Seed>& seeds, const std::vector<Eigen::MatrixXd>& V,
                   const std::vector<Eigen::MatrixXi>& F, std::array<Eigen::MatrixXd, 4>& SV,
                   std::array<Eigen::MatrixXi, 4>& SF)
{
    const auto k = patches(BV, seeds);
    std::array<std::vector<Eigen::MatrixXd>, 4> PV;
    std::array<std::vector<Eigen::MatrixXi>, 4> PF;
    for (size_t i = 0; i < k.size(); ++i) {
        PV[k[i]].push_back(V[i]);
        PF[k[i]].push_back(F[i]);
    }
    for (int j = 0; j < 4; ++j)
        if (!PV[j].empty()) remesh(PV[j], PF[j], 0, SV[j], SF[j]);
}

} // namespace pavel
