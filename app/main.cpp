#include "deform.hpp"
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
                       "without it, the first element is seeded isotropically")
            ->check(CLI::ExistingFile);
    app.add_option("--coverage", iso.coverage, "Isotropic seeding: footprint area of all placed elements / base area")
        ->check(CLI::PositiveNumber)
        ->excludes(seeds_opt)
        ->capture_default_str();
    app.add_option("--lloyd", iso.lloyd, "Isotropic seeding: CVT Lloyd iterations")
        ->check(CLI::NonNegativeNumber)
        ->excludes(seeds_opt)
        ->capture_default_str();
    app.add_option("--newton", iso.newton, "Isotropic seeding: CVT Newton iterations")
        ->check(CLI::NonNegativeNumber)
        ->excludes(seeds_opt)
        ->capture_default_str();
    app.add_option("--voxel-size", h, "Voxel edge length")->required()->check(CLI::PositiveNumber);
    app.add_option("-a", p.a, "Recovery function: peak at a·d")->capture_default_str();
    app.add_option("-b", p.b, "Recovery function: zero beyond b·d")->capture_default_str();
    app.add_option("--order", p.order, "Fast marching order")->check(CLI::Range(1, 2))->capture_default_str();
    app.add_option("-o,--output", out, "Output directory")->capture_default_str();
    CLI11_PARSE(app, argc, argv);

    Eigen::MatrixXd BV;
    Eigen::MatrixXi BF;
    read_mesh(base_path, BV, BF);
    std::vector<Eigen::MatrixXd> EV(element_paths.size());
    std::vector<Eigen::MatrixXi> EF(element_paths.size());
    for (size_t k = 0; k < element_paths.size(); ++k) read_mesh(element_paths[k], EV[k], EF[k]);

    std::ifstream in(seeds_path);
    const auto seeds = seeds_path.empty()
                           ? pavel::isotropic_seeds(BV, BF, pavel::footprint_area(EV[0], EF[0], h), iso)
                           : pavel::read_seeds(in);
    std::cout << seeds.size() << " seeds\n";

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
    std::filesystem::create_directories(out);
    igl::write_triangle_mesh((out / "elements.obj").string(), VC, FC);
}
