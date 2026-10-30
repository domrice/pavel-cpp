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
m=assets/meshes
$pavel --base $m/hemisphere.obj --elements $m/ball.obj --voxel-size 0.005                    # isotropic seeding, coverage 1.4
$pavel --base $m/hemisphere.obj --elements $m/ball.obj --voxel-size 0.005 --offset --shells  # offset-surface seeding, printable shells
$pavel --base base.obj --elements a.obj b.obj --voxel-size 0.05 \
       --stripes 1.5 1.5 90 --rotation alternate --element-choice alternate                  # stripe seeding
```

`assets/meshes` holds two test bases, `hemisphere.obj` and `saddle.obj` (radius 1, +z up), and one element, `ball.obj` (radius 0.1, sunk 0.01 below z = 0). Both bases are open surfaces. The inside test still clips the elements against them, but the pipeline is meant for closed bases. Model each element in its own frame: +z is its symmetry axis, and z = 0 is the footprint plane that sits on the base. Whatever lies below z = 0 ends up inside the base and is lost for good. A seeds file has one line per element, `x y z nx ny nz angle element scale`. `--voxel-size` sets the resolution: the paper's per-element grids have diagonals of roughly 120–200 voxels. For shells, orient the base so that its symmetry axis is +z. `$pavel --help` lists every option. Options can also come from a TOML/INI file passed with `--config`.

Results go to `out/` (`-o` to change it): `seeds.txt` (the seeds used, in the `--seeds` format), `elements.obj` (the deformed elements), `remeshed.obj` (elements + base, CVT-remeshed), and with `--shells` also `shell0..3.obj`.

## Placing seeds by hand

Configure with `-DPAVEL_VIEWER=ON`, build, and add `--place`. An empty seeds file starts from scratch. Without `--seeds`, the window opens on the automatic seeds (isotropic or `--stripes`).

```bash
touch empty.txt
$pavel --base $m/hemisphere.obj --elements $m/ball.obj --voxel-size 0.005 --seeds empty.txt --place
```

In the polyscope window, drag to rotate, right-drag to pan, and scroll to zoom. A click on the base adds a seed there, carrying the previous seed's element, angle and scale. The panel edits only the newest seed: `element` (index into `--elements`), `angle` (degrees about the normal) and `scale`. Undo removes the newest seed. Deform previews the deformed elements on the current seeds, and the window freezes until it is done. Any edit after that drops the preview. Closing the window runs the full pipeline on the edited seeds and writes them to `out/seeds.txt`. To continue editing later, pass that file back with `--seeds`.

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
assets/     test meshes: two bases and an element
cmake/      pinned third-party dependencies
out/        run outputs
```
