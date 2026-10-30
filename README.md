# PAVEL C++

[![License](https://img.shields.io/github/license/domrice/pavel-cpp?label=license)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599c?logo=cplusplus&logoColor=white)](https://en.cppreference.com/w/cpp/20)
[![CMake 3.25+](https://img.shields.io/badge/CMake-3.25%2B-064f8c?logo=cmake&logoColor=white)](https://cmake.org/)
[![Paper-DOI](https://img.shields.io/badge/paper-10.1145%2F3502802-215CAF)](https://doi.org/10.1145/3502802)

A C++20 reimplementation of PAVEL (Fanni, Pellacini, Scateni, Giachetti, ACM TOG 41(2), 2022). It places volumetric decorative elements on a base mesh and deforms them, preserving their volume, until they pack against each other like clay. Output is a smooth mesh or a set of printable shells. The authors' code is not public, so this one follows the paper as literally as possible. Ambiguous aspects of the paper are interpreted.

## Quick start

Needs CMake ≥ 3.25, Ninja, a C++20 compiler, and git. The first configure fetches libigl + Eigen, geometry-central, Geogram, CLI11 and Catch2, and polyscope too when the viewer is on. Nothing is installed system-wide. Needs [uv](https://docs.astral.sh/uv/getting-started/installation/) only to regenerate the test fixtures.

```bash
cmake --preset release                       # add -DPAVEL_VIEWER=ON for the placement viewer
cmake --build --preset release
ctest --preset release

pavel=./build/release/app/pavel
$pavel --base base.obj --elements ball.obj --voxel-size 0.05                         # isotropic seeding, coverage 1.4
$pavel --base base.obj --elements ball.obj --voxel-size 0.05 --offset --shells       # offset-surface seeding, printable shells
$pavel --base base.obj --elements a.obj b.obj --voxel-size 0.05 \
       --stripes 1.5 1.5 90 --rotation alternate --element-choice alternate          # stripe seeding
$pavel --base base.obj --elements ball.obj --voxel-size 0.05 --seeds seeds.txt --place   # edit seeds in the viewer
```

No meshes ship with the repo. Model each element in its own frame: +z is its symmetry axis, and z = 0 is the footprint plane that sits on the base. A seeds file has one line per element, `x y z nx ny nz angle element scale`. `--voxel-size` sets the resolution: the paper's per-element grids have diagonals of roughly 120–200 voxels. For shells, orient the base so that its symmetry axis is +z. `$pavel --help` lists every option. Options can also come from a TOML/INI file passed with `--config`.

Results go to `out/` (`-o` to change it): `seeds.txt` (the seeds used, in the `--seeds` format), `elements.obj` (the deformed elements), `remeshed.obj` (elements + base, CVT-remeshed), and with `--shells` also `shell0..3.obj`.

Test fixtures:

```bash
uv run tools/fmm_reference.py   # scikit-fmm arrival times → tests/data/fmm
```

## Layout

```
src/        header-only library pavel::core, one header per pipeline stage: seeding, voxel_grid,
            voxelize, fast_marching, recovery_function, deform, output, mesh_convert
app/        CLI `pavel` (main.cpp), polyscope placement viewer (viewer.hpp)
tests/      Catch2 tests, one file per module; scikit-fmm reference data in tests/data/
tools/      uv single-file scripts (PEP 723 + lock file)
cmake/      pinned third-party dependencies
docs/       the paper (TOG text, arXiv TeX) and the roadmap
out/        run outputs
```
