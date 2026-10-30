#pragma once

#include "mesh_convert.hpp"
#include "voxelize.hpp"

#include <Eigen/Geometry>
#include <geogram/voronoi/CVT.h>
#include <geometrycentral/surface/direction_fields.h>
#include <geometrycentral/surface/stripe_patterns.h>
#include <geometrycentral/surface/surface_mesh_factories.h>
#include <igl/barycentric_coordinates.h>
#include <igl/doublearea.h>
#include <igl/per_vertex_normals.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/signed_distance.h>

#include <cmath>
#include <istream>
#include <numbers>
#include <ostream>
#include <queue>
#include <random>
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

// Seeds in the format read_seeds reads.
inline void write_seeds(std::ostream& out, const std::vector<Seed>& seeds)
{
    for (const auto& s : seeds)
        out << s.position.transpose() << ' ' << s.normal.transpose() << ' ' << s.angle << ' ' << s.element << ' '
            << s.scale << '\n';
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

// The angle (degrees) at which place() turns the element's +x along the tangent t at normal n, read off the unturned
// placement of +x and +y.
inline double angle_toward(const Eigen::Vector3d& n, const Eigen::Vector3d& t)
{
    const Eigen::Vector2d c = place(Eigen::MatrixXd{{1, 0, 0}, {0, 1, 0}}, {Eigen::Vector3d::Zero(), n, 0, 0, 1}) * t;
    return std::atan2(c.y(), c.x()) * 180 / std::numbers::pi;
}

// PAPER: §4, the footprint of an element (V, F in its frame) is its z-section of maximal area, measured on its voxels
// of size h. Its height (§5) is the offset d of offset-surface seeding. Automatic seeding chooses among the elements
// V[k], F[k] uniformly (§5), so it returns their mean.
struct Footprint {
    double area, height;
};

inline Footprint footprint(const std::vector<Eigen::MatrixXd>& V, const std::vector<Eigen::MatrixXi>& F, double h)
{
    Footprint fp{0, 0};
    for (size_t k = 0; k < V.size(); ++k) {
        const auto g = voxelize(V[k], F[k], h, element_box(V[k]));
        const auto layers = g.inside.cast<int>().reshaped(Eigen::Index(g.size.x()) * g.size.y(), g.size.z());
        Eigen::Index z;
        const int n = layers.colwise().sum().maxCoeff(&z);
        fp.area += h * h * n / V.size();
        fp.height += h * (g.offset.z() + z) / V.size();
    }
    return fp;
}

// PAPER: §5, the isosurface OV, OF at signed distance d from the base V, F: marching cubes of the signed distance
// (sign by fast winding number) sampled at the lattice voxel centers of size h.
// ponytail: the distance is computed over the base's whole bounding box; bound it to a band around d (igl's
// lower/upper_bound overload) if this dominates on large bases.
inline void offset_surface(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double d, double h, Eigen::MatrixXd& OV,
                           Eigen::MatrixXi& OF)
{
    const Eigen::Vector3d margin = Eigen::Vector3d::Constant(std::abs(d) + 2 * h);
    const auto g = lattice_grid(h, {V.colwise().minCoeff().transpose() - margin,
                                    V.colwise().maxCoeff().transpose() + margin});
    const Eigen::MatrixXd P = g.centers();
    Eigen::VectorXd S;
    Eigen::VectorXi I;
    Eigen::MatrixXd C, N;
    igl::signed_distance(P, V, F, igl::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER, S, I, C, N);
    igl::marching_cubes(S, P, g.size.x(), g.size.y(), g.size.z(), d, OV, OF);
}

// Seeds at the points P projected onto the base V, F (closest point), each with the base's vertex normals
// interpolated there, placing element 0 unturned and unscaled.
inline std::vector<Seed> surface_seeds(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& P)
{
    Eigen::MatrixXd C, B, N;
    Eigen::VectorXd D;
    Eigen::VectorXi I;
    igl::point_mesh_squared_distance(P, V, F, D, I, C);
    using Eigen::placeholders::all;
    igl::barycentric_coordinates(C, V(F(I, 0), all), V(F(I, 1), all), V(F(I, 2), all), B);
    igl::per_vertex_normals(V, F, N);
    std::vector<Seed> seeds;
    for (Eigen::Index i = 0; i < P.rows(); ++i) {
        const Eigen::Vector3d normal = B.row(i) * N(F.row(I[i]), all);
        seeds.push_back({C.row(i), normal.normalized(), 0, 0, 1});
    }
    return seeds;
}

// PAPER: §4, coverage about 1.4 (1.2 barely deforms, 1.6 distorts). The iteration counts are GEO::remesh_smooth's.
struct IsotropicParams {
    double coverage = 1.4;
    int lloyd = 5, newton = 30;
};

// Isotropic seeding (§5) on the base V, F: coverage = N · footprint / base area fixes the number N of seeds, which
// Geogram's CVT (Yan et al. 2009) spreads over the sampling surface SV, SF (the base, or its offset surface) with
// lloyd then newton iterations; fewer give a less regular, "hand-made" look. Each seed is projected onto the base
// (surface_seeds; for the offset surface, the closest point is the projection along its normal).
inline std::vector<Seed> isotropic_seeds(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& SV,
                                         const Eigen::MatrixXi& SF, double footprint, const IsotropicParams& p)
{
    Eigen::VectorXd A;
    igl::doublearea(V, F, A);
    const auto n = GEO::index_t(std::lround(p.coverage * A.sum() / 2 / footprint));
    GEO::Mesh M;
    to_geogram(SV, SF, M);
    GEO::CentroidalVoronoiTesselation cvt(&M);
    cvt.compute_initial_sampling(n);
    cvt.Lloyd_iterations(p.lloyd);
    if (p.newton > 0) cvt.Newton_iterations(p.newton); // as in GEO::remesh_smooth

    return surface_seeds(V, F,
                         Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>>(
                             cvt.embedding(0), n, 3));
}

// PAPER: §5, the stripe seeds' connectivity is "estimated": seed j neighbours seed i along a stripe family when it lies
// within 1.5 gaps (the distance between consecutive crossings along that family) of i, at most 20° off the family's
// direction at i (T1, T2). A breadth-first walk over these links flips the main directions T1 to the same sign across
// every link, and returns each seed's index along its main stripe mod 2: it flips along the main stripes and is kept
// along the second ones.
// ponytail: O(n²) neighbour search, and links that disagree (around branch points) lose to the first visit; a k-d tree
// once seed counts reach the 10⁴s.
inline std::vector<int> orient_stripes(const std::vector<Seed>& seeds, std::vector<Eigen::Vector3d>& T1,
                                       const std::vector<Eigen::Vector3d>& T2, double gap1, double gap2)
{
    const double c = std::cos(20 * std::numbers::pi / 180);
    std::vector<int> parity(seeds.size(), -1);
    std::queue<size_t> q;
    for (size_t r = 0; r < seeds.size(); ++r) {
        if (parity[r] != -1) continue;
        parity[r] = 0;
        q.push(r);
        while (!q.empty()) {
            const size_t i = q.front();
            q.pop();
            for (size_t j = 0; j < seeds.size(); ++j) {
                if (parity[j] != -1) continue;
                const Eigen::Vector3d e = seeds[j].position - seeds[i].position;
                const bool along1 = e.norm() < 1.5 * gap1 && std::abs(e.normalized().dot(T1[i])) > c;
                const bool along2 = e.norm() < 1.5 * gap2 && std::abs(e.normalized().dot(T2[i])) > c;
                if (!along1 && !along2) continue;
                if (T1[j].dot(T1[i]) < 0) T1[j] = -T1[j];
                parity[j] = parity[i] ^ along1;
                q.push(j);
            }
        }
    }
    return parity;
}

// Stripe seeding (§5) on the base V, F: two Knöppel et al. (2015) stripe patterns, the main stripes distance apart
// along the smoothest line field, the second ones crossing them at angle (degrees) and spaced so that consecutive
// crossings along a main stripe lie spacing apart. Seeds sit at the crossings (none in the patterns' singular faces).
// Their normal interpolates the base's vertex normals, and they place element 0 unscaled, turned so that its +x runs
// along the main stripe. parity[i] is seed i's index along its main stripe mod 2 (orient_stripes).
// ponytail: a crossing on a shared edge may be found in both faces or neither (round-off); dedupe if it shows.
inline std::vector<Seed> stripe_seeds(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, double distance,
                                      double spacing, double angle, std::vector<int>& parity)
{
    namespace gcs = geometrycentral::surface;
    const double alpha = angle * std::numbers::pi / 180;
    const auto [mesh, geom] = gcs::makeManifoldSurfaceMeshAndGeometry(V, F);
    const auto X1 = gcs::computeSmoothestVertexDirectionField(*geom, 2);
    gcs::VertexData<geometrycentral::Vector2> X2(*mesh);
    for (const auto v : mesh->vertices()) X2[v] = X1[v] * geometrycentral::Vector2::fromAngle(2 * alpha); // 2-RoSy
    // Frequencies are stripes per unit length; the second stripes lie spacing·sin(alpha) apart.
    const auto [t1, s1, b1] = gcs::computeStripePattern(*geom, gcs::VertexData<double>(*mesh, 1 / distance), X1);
    const auto [t2, s2, b2] =
        gcs::computeStripePattern(*geom, gcs::VertexData<double>(*mesh, 1 / (spacing * std::sin(alpha))), X2);

    Eigen::MatrixXd N;
    igl::per_vertex_normals(V, F, N);
    std::vector<Seed> seeds;
    std::vector<Eigen::Vector3d> T1, T2;
    for (const auto f : mesh->faces()) {
        if (s1[f] != 0 || s2[f] != 0 || b1[f] != 0) continue; // b2 == b1: X2 is X1 turned by a constant
        const auto h0 = f.halfedge(), h1 = h0.next(), h2 = h1.next();
        const Eigen::Vector3i v(h0.vertex().getIndex(), h1.vertex().getIndex(), h2.vertex().getIndex());
        // Both phases, in turns, are affine over the face: a = a0 + l0 (a1 - a0) + l1 (a2 - a0), likewise b.
        const Eigen::Vector3d a = Eigen::Vector3d(t1[h0.corner()], t1[h1.corner()], t1[h2.corner()]) / (2 * std::numbers::pi);
        const Eigen::Vector3d b = Eigen::Vector3d(t2[h0.corner()], t2[h1.corner()], t2[h2.corner()]) / (2 * std::numbers::pi);
        const Eigen::Matrix2d M{{a[1] - a[0], a[2] - a[0]}, {b[1] - b[0], b[2] - b[0]}};
        const Eigen::Vector3d p0 = V.row(v[0]), e1 = V.row(v[1]) - V.row(v[0]), e2 = V.row(v[2]) - V.row(v[0]);
        for (int k = int(std::ceil(a.minCoeff())); k <= a.maxCoeff(); ++k)
            for (int m = int(std::ceil(b.minCoeff())); m <= b.maxCoeff(); ++m) {
                const Eigen::Vector2d l = M.inverse() * Eigen::Vector2d(k - a[0], m - b[0]);
                if (!(l.minCoeff() >= 0 && l.sum() <= 1)) continue; // also rejects NaN from a degenerate M
                const Eigen::Vector3d n = Eigen::RowVector3d(1 - l.sum(), l[0], l[1]) * N(v, Eigen::placeholders::all);
                seeds.push_back({p0 + l[0] * e1 + l[1] * e2, n.normalized(), 0, 0, 1});
                // Each family's isoline direction: where its phase stays constant.
                T1.push_back(((a[2] - a[0]) * e1 - (a[1] - a[0]) * e2).normalized());
                T2.push_back(((b[2] - b[0]) * e1 - (b[1] - b[0]) * e2).normalized());
            }
    }

    parity = orient_stripes(seeds, T1, T2, spacing, distance / std::sin(alpha));
    for (size_t i = 0; i < seeds.size(); ++i) seeds[i].angle = angle_toward(seeds[i].normal, T1[i]);
    return seeds;
}

// §5, see PerturbParams.
enum class Variation { none, random, alternate };

// PAPER: §5, element perturbation. Rotation: random adds a uniform turn about the normal, alternate adds 180° every
// other seed along the main stripe (none keeps the seeding's orientation: unturned, or along the main stripe). Element:
// random picks one of the loaded elements uniformly, alternate takes parity mod their count. Scale: a uniform factor
// in 1 ± scale_jitter, which the paper caps at 0.1.
struct PerturbParams {
    Variation rotation = Variation::none, element = Variation::random;
    double scale_jitter = 0;
    unsigned random_seed = 0;
};

// Perturbs automatically placed seeds, choosing among the first `elements` loaded elements; parity comes from
// stripe_seeds and is read only to alternate.
inline void perturb(std::vector<Seed>& seeds, const std::vector<int>& parity, int elements, const PerturbParams& p)
{
    std::mt19937 rng(p.random_seed);
    std::uniform_real_distribution<double> turn(0, 360), jitter(-p.scale_jitter, p.scale_jitter);
    std::uniform_int_distribution<int> pick(0, elements - 1);
    for (size_t i = 0; i < seeds.size(); ++i) {
        Seed& s = seeds[i];
        if (p.rotation == Variation::random) s.angle += turn(rng);
        if (p.rotation == Variation::alternate) s.angle += 180 * parity[i];
        if (p.element == Variation::random) s.element = pick(rng);
        if (p.element == Variation::alternate) s.element = parity[i] % elements;
        s.scale *= 1 + jitter(rng);
    }
}

} // namespace pavel
