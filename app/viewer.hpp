#pragma once

#include "output.hpp"

#include <imgui.h>
#include <polyscope/pick.h>
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>

// Shows the meshes V[k], F[k] as one structure, or removes it when there are none.
static void show(const std::string& name, const std::vector<Eigen::MatrixXd>& V, const std::vector<Eigen::MatrixXi>& F)
{
    if (V.empty()) return polyscope::removeStructure(name);
    Eigen::MatrixXd CV;
    Eigen::MatrixXi CF;
    igl::combine(V, F, CV, CF);
    polyscope::registerSurfaceMesh(name, CV, CF);
}

// Shows the elements EV, EF placed at the seeds, and drops the deformed preview, which they outdate.
static void show_seeds(const std::vector<Eigen::MatrixXd>& EV, const std::vector<Eigen::MatrixXi>& EF,
                       const std::vector<pavel::Seed>& seeds)
{
    std::vector<Eigen::MatrixXd> V;
    std::vector<Eigen::MatrixXi> F;
    for (const auto& s : seeds) {
        V.push_back(pavel::place(EV.at(s.element), s));
        F.push_back(EF.at(s.element));
    }
    show("seeds", V, F);
    polyscope::removeStructure("deformed");
}

// Manual placement (§5) in a polyscope window over the base BV, BF: a click on the base adds a seed there, the panel
// sets the last seed's element, angle and scale, a new seed starts from the last one's, Undo removes it, and Deform
// previews the pipeline's deformed elements. seeds holds the edited seeds when the window closes.
static void edit_seeds(const Eigen::MatrixXd& BV, const Eigen::MatrixXi& BF, const std::vector<Eigen::MatrixXd>& EV,
                       const std::vector<Eigen::MatrixXi>& EF, double h, const pavel::DeformParams& p,
                       std::vector<pavel::Seed>& seeds)
{
    polyscope::init();
    const auto base = polyscope::registerSurfaceMesh("base", BV, BF);
    show_seeds(EV, EF, seeds);
    polyscope::state::userCallback = [&] {
        const auto pick = polyscope::getSelection();
        bool changed = pick.isHit && pick.structure == base;
        if (changed) {
            auto s = pavel::surface_seeds(BV, BF, Eigen::RowVector3d(pick.position.x, pick.position.y,
                                                                      pick.position.z))[0];
            if (!seeds.empty()) s = {s.position, s.normal, seeds.back().angle, seeds.back().element, seeds.back().scale};
            seeds.push_back(s);
            polyscope::resetSelection();
        }
        ImGui::TextUnformatted("Click the base to add a seed.");
        if (!seeds.empty()) {
            auto& s = seeds.back();
            changed |= ImGui::SliderInt("element", &s.element, 0, int(EV.size()) - 1);
            changed |= ImGui::InputDouble("angle", &s.angle, 15);
            changed |= ImGui::InputDouble("scale", &s.scale, 0.05);
        }
        if (ImGui::Button("Undo") && !seeds.empty()) {
            seeds.pop_back();
            changed = true;
        }
        if (changed) show_seeds(EV, EF, seeds);
        ImGui::SameLine();
        if (ImGui::Button("Deform") && !seeds.empty()) {
            std::vector<Eigen::MatrixXd> V;
            std::vector<Eigen::MatrixXi> F;
            pavel::deformed_meshes(BV, BF, EV, EF, seeds, h, p, V, F);
            show("deformed", V, F);
            polyscope::getSurfaceMesh("seeds")->setEnabled(false);
        }
        ImGui::Text("%zu seeds", seeds.size());
    };
    polyscope::show();
}
