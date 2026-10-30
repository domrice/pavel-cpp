# PAVEL C++

[![License](https://img.shields.io/github/license/domrice/pavel-cpp?label=license)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599c?logo=cplusplus&logoColor=white)](https://en.cppreference.com/w/cpp/20)
[![CMake 3.25+](https://img.shields.io/badge/CMake-3.25%2B-064f8c?logo=cmake&logoColor=white)](https://cmake.org/)
[![Paper-DOI](https://img.shields.io/badge/paper-10.1145%2F3502802-215CAF)](https://doi.org/10.1145/3502802)

A C++20 reimplementation of **PAVEL: Decorative Patterns with Packed Volumetric Elements** by Fanni et al. It places volumetric elements on a base mesh and deforms them, preserving their volume, until they pack against each other like clay. The authors' code is not public, so this one follows the paper as literally as possible; where the paper is ambiguous, [AGENTS.md](AGENTS.md) lists the interpretation taken.

## Build

Needs CMake ≥ 3.25, Ninja, a C++20 compiler and git. The first configure fetches all dependencies into `build/`; nothing is installed system-wide.

```bash
cmake --preset release            # configure, download and build the dependencies (first run only)
cmake --build --preset release    # build the pavel CLI and the tests
ctest --preset release            # run the tests (optional)
```

Add `-DPAVEL_VIEWER=ON` to the first line for manual placement (`--place`).

## Run on your own meshes

```bash
pavel=./build/release/app/pavel
$pavel --base base.obj --elements a.obj b.obj --voxel-size 0.01
```

Inputs:

- `--base`: a closed triangle mesh (OBJ, OFF, STL or PLY). For `--shells`, orient its symmetry axis along +z.
- `--elements`: one mesh per element, in the base's units, with +z as its symmetry axis and z = 0 as the plane that sits on the base. Anything below z = 0 is lost inside the base.
- `--voxel-size`: in mesh units. Each element's grid is 1.5× its bounding box; aim for a grid diagonal of 120–200 voxels, as in the paper. Start coarse.

Seeding, one of:

- default: isotropic, `--coverage` sets the density (footprint area of all elements / base area, default 1.4), `--offset` seeds on the offset surface.
- `--stripes D L α`: stripes D apart, seeds L apart along them, crossing angle α in degrees.
- `--seeds f.txt`: one line per element, `x y z nx ny nz angle element scale` (angle in degrees about the normal, element as an index into `--elements`). Used as is, without perturbation.

Perturbation, for isotropic and stripe seeds:

- `--rotation`: turn about the normal, `none` (default), `random`, or `alternate` (180° every other seed, `--stripes` only).
- `--element-choice`: which of `--elements`, `none` (the first), `random` (default), or `alternate` (`--stripes` only).
- `--scale-jitter`: uniform scale factor in 1 ± this, from 0 (default) to 0.1.
- `--random-seed`: seed of the random numbers.

Outputs, in `out/` (`-o` to change it):

- `seeds.txt`: the seeds used, reusable with `--seeds`.
- `elements.obj`: the deformed elements.
- `remeshed.obj`: elements + base, CVT-remeshed.
- `shell0..3.obj`: printable quadrants, with `--shells`.

`$pavel --help` lists every option; `--config f.toml` reads them from a file.

Examples on the bundled meshes (`assets/meshes`: open bases `hemisphere.obj`, `saddle.obj` of radius 1; element `ball.obj` of radius 0.1):

```bash
m=assets/meshes
$pavel --base $m/hemisphere.obj --elements $m/ball.obj --voxel-size 0.005
$pavel --base $m/hemisphere.obj --elements $m/ball.obj --voxel-size 0.005 --offset --shells
$pavel --base $m/saddle.obj --elements $m/ball.obj --voxel-size 0.005 --stripes 0.2 0.2 90 --rotation alternate
```

## Placing seeds by hand

Build with `-DPAVEL_VIEWER=ON` and add `--place`. Without `--seeds` the window opens on the automatic seeds; an empty seeds file starts from none.

```bash
touch empty.txt
$pavel --base $m/hemisphere.obj --elements $m/ball.obj --voxel-size 0.005 --seeds empty.txt --place
```

Click the base to add a seed (it copies the previous seed's element, angle and scale). The panel edits the newest seed; Undo removes it; Deform previews the result. Closing the window runs the full pipeline and writes `out/seeds.txt`, which you can pass back with `--seeds` to keep editing.

## Citation

I did not contribute to this paper, but please cite the authors' work when you use this code.

```
@article{fanniPAVELDecorativePatterns2022,
  title = {{{PAVEL}}: {{Decorative Patterns}} with {{Packed Volumetric Elements}}},
  author = {Fanni, Filippo Andrea and Pellacini, Fabio and Scateni, Riccardo and Giachetti, Andrea},
  year = 2022,
  month = apr,
  journal = {ACM Transactions on Graphics},
  volume = {41},
  number = {2},
  pages = {1--15},
  publisher = {Association for Computing Machinery (ACM)},
  doi = {10.1145/3502802}
}
```
