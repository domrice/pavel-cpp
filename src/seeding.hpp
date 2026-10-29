#pragma once

#include <Eigen/Geometry>

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

} // namespace pavel
