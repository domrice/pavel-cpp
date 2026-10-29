#include "deform.hpp"
#include "output.hpp"
#include "seeding.hpp"
#include "voxelize.hpp"

#include <CLI/CLI.hpp>
#include <igl/combine.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>

#include <filesystem>
#include <fstream>
#include <iostream>

static void read_mesh(const std::string& path, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
    if (!igl::read_triangle_mesh(path, V, F)) throw std::runtime_error("cannot read " + path);
}

int main(int argc, char** argv)
{
    CLI::App app{"PAVEL: decorative patterns with packed volumetric elements"};
    argv = app.ensure_utf8(argv);
    app.set_config("--config");
    std::string base_path, seeds_path;
    std::vector<std::string> element_paths;
    double h;
    pavel::IsotropicParams iso;
    pavel::DeformParams p;
    std::filesystem::path out = "out";
    app.add_option("--base", base_path, "Base mesh")->required()->check(CLI::ExistingFile);
    app.add_option("--elements", element_paths, "Element meshes: +z is the symmetry axis, z = 0 the footprint plane")
        ->required()
        ->check(CLI::ExistingFile);
    const auto seeds_opt =
        app.add_option("--seeds", seeds_path,
                       "Seeds, one per line: x y z nx ny nz angle(degrees) element(index into --elements) scale; "
                       "without it (and --stripes), the elements are seeded isotropically")
            ->check(CLI::ExistingFile);
    std::vector<double> stripes;
    const auto stripes_opt =
        app.add_option("--stripes", stripes,
                       "Stripe seeding: distance between main stripes, seed spacing along them, crossing angle "
                       "(degrees)")
            ->expected(3)
            ->check(CLI::PositiveNumber)
            ->excludes(seeds_opt);
    const auto iso_group = app.add_option_group("Isotropic seeding");
    iso_group->excludes(seeds_opt);
    iso_group->excludes(stripes_opt);
    bool offset = false;
    iso_group->add_flag("--offset", offset, "Seed on the offset surface at the footprint's height");
    iso_group->add_option("--coverage", iso.coverage, "Footprint area of all placed elements / base area")
        ->check(CLI::PositiveNumber)
        ->capture_default_str();
    iso_group->add_option("--lloyd", iso.lloyd, "CVT Lloyd iterations")
        ->check(CLI::NonNegativeNumber)
        ->capture_default_str();
    iso_group->add_option("--newton", iso.newton, "CVT Newton iterations")
        ->check(CLI::NonNegativeNumber)
        ->capture_default_str();
    const auto perturb_group = app.add_option_group("Perturbation");
    perturb_group->excludes(seeds_opt);
    pavel::PerturbParams pp;
    const std::map<std::string, pavel::Variation> variation{
        {"none", pavel::Variation::none}, {"random", pavel::Variation::random}, {"alternate", pavel::Variation::alternate}};
    perturb_group->add_option("--rotation", pp.rotation,
                   "Turn about the normal: none, random, or alternate (180° every other seed along the main stripe, "
                   "--stripes only); default none")
        ->transform(CLI::CheckedTransformer(variation));
    perturb_group->add_option("--element-choice", pp.element,
                   "Element among --elements: none (the first), random, or alternate (along the main stripe, "
                   "--stripes only); default random")
        ->transform(CLI::CheckedTransformer(variation));
    perturb_group->add_option("--scale-jitter", pp.scale_jitter, "Uniform scale factor in 1 ± this")
        ->check(CLI::Range(0.0, 0.1))
        ->capture_default_str();
    perturb_group->add_option("--random-seed", pp.random_seed, "Seed of the perturbations' random numbers")
        ->capture_default_str();
    app.add_option("--voxel-size", h, "Voxel edge length")->required()->check(CLI::PositiveNumber);
    app.add_option("-a", p.a, "Recovery function: peak at a·d")->capture_default_str();
    app.add_option("-b", p.b, "Recovery function: zero beyond b·d")->capture_default_str();
    app.add_option("--order", p.order, "Fast marching order")->check(CLI::Range(1, 2))->capture_default_str();
    unsigned remesh_points = 0;
    app.add_option("--remesh-points", remesh_points,
                   "CVT remesh of elements + base: number of points; 0 keeps the merged mesh's vertex count")
        ->capture_default_str();
    app.add_option("-o,--output", out, "Output directory")->capture_default_str();
    CLI11_PARSE(app, argc, argv);
    if (stripes.empty() && (pp.rotation == pavel::Variation::alternate || pp.element == pavel::Variation::alternate))
        return app.exit(CLI::ValidationError("alternate", "needs --stripes"));

    Eigen::MatrixXd BV;
    Eigen::MatrixXi BF;
    read_mesh(base_path, BV, BF);
    std::vector<Eigen::MatrixXd> EV(element_paths.size());
    std::vector<Eigen::MatrixXi> EF(element_paths.size());
    for (size_t k = 0; k < element_paths.size(); ++k) read_mesh(element_paths[k], EV[k], EF[k]);

    std::vector<pavel::Seed> seeds;
    std::vector<int> parity;
    if (!seeds_path.empty()) {
        std::ifstream in(seeds_path);
        seeds = pavel::read_seeds(in);
    } else {
        if (!stripes.empty()) {
            seeds = pavel::stripe_seeds(BV, BF, stripes[0], stripes[1], stripes[2], parity);
        } else {
            const auto fp = pavel::footprint(EV, EF, h);
            Eigen::MatrixXd SV = BV;
            Eigen::MatrixXi SF = BF;
            if (offset) pavel::offset_surface(BV, BF, fp.height, h, SV, SF);
            seeds = pavel::isotropic_seeds(BV, BF, SV, SF, fp.area, iso);
        }
        pavel::perturb(seeds, parity, int(EV.size()), pp);
    }
    std::cout << seeds.size() << " seeds\n";
    std::filesystem::create_directories(out);
    std::ofstream seeds_out(out / "seeds.txt");
    pavel::write_seeds(seeds_out, seeds);

    // ponytail: the base's winding-number tree is rebuilt per element; precompute it once if that shows up in a
    // profile.
    std::vector<pavel::VoxelGrid> elements;
    std::vector<pavel::Mask> base;
    for (const auto& s : seeds) {
        const Eigen::MatrixXd V = pavel::place(EV.at(s.element), s);
        const auto box = pavel::element_box(V);
        elements.push_back(pavel::voxelize(V, EF.at(s.element), h, box));
        base.push_back(pavel::voxelize(BV, BF, h, box).inside);
    }
    const auto deformed = pavel::deform(elements, base, p);

    std::vector<Eigen::MatrixXd> V(deformed.size());
    std::vector<Eigen::MatrixXi> F(deformed.size());
    for (size_t i = 0; i < deformed.size(); ++i) {
        if (deformed[i].inside.count() < elements[i].inside.count())
            std::cout << "element " << i << " recovered " << deformed[i].inside.count() << " of "
                      << elements[i].inside.count() << " voxels\n";
        pavel::marching_cubes(deformed[i], V[i], F[i]);
    }
    Eigen::MatrixXd VC;
    Eigen::MatrixXi FC;
    igl::combine(V, F, VC, FC);
    igl::write_triangle_mesh((out / "elements.obj").string(), VC, FC);

    pavel::remesh({VC, BV}, {FC, BF}, remesh_points, VC, FC);
    igl::write_triangle_mesh((out / "remeshed.obj").string(), VC, FC);
}
