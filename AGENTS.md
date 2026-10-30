# AGENTS.md

Guidance for coding agents (and people) working on this repository. Usage is in [README.md](README.md).

## Project

A C++20 reimplementation of **PAVEL** (Fanni, Pellacini, Scateni, Giachetti, *PAVEL: Decorative Patterns with Packed Volumetric Elements*, ACM TOG 41(2), 2022, doi:10.1145/3502802; preprint arXiv:2102.01029). It aims to be as faithful as possible to the authors' implementation, which is not public. Per the paper, they used Python + scikit-fmm (fast marching) + trimesh (IO, voxelization) + Geogram (CVT seeding and CVT remeshing).

The paper is the spec. The final TOG version is authoritative and is the only one with §8 and with phase A/B of the deformation. The arXiv TeX has the readable equations (Eqs. 1 and 2 are garbled in the TOG text extraction). When the paper is ambiguous, implement the most literal reading, mark it with `// PAPER:` + section, and add it to "Open interpretations" below. Do not silently "improve" the algorithm.

## Commands

```sh
cmake --preset release                 # configure (deps fetched on first run); add -DPAVEL_VIEWER=ON for --place
cmake --build --preset release         # build
ctest --preset release                 # all tests
./build/release/tests/pavel_tests "[recovery]"   # one Catch2 tag
./build/release/tests/pavel_tests "[visual]"     # hidden: writes two-ball meshes to out/
./build/release/app/pavel --help       # CLI
uv run tools/fmm_reference.py          # regenerate tests/data/fmm from scikit-fmm (locked in tools/*.lock)
```

Presets: `debug` (ASan/UBSan, deps included) and `release` (`-O3`). Time things in `release` only: the deformation handles about 10⁸ voxels. Each preset fetches and builds its own deps on first configure.

Known harmless noise: `ld: warning: search path …/geogram-src/build/release/lib not found`.

## Layout

```
cmake/Dependencies.cmake   third-party deps via FetchContent, pinned tags (Geogram by tarball SHA256)
src/                       header-only INTERFACE library pavel::core, one header per pipeline stage, namespace pavel
app/                       CLI `pavel` (main.cpp, CLI11, --config for TOML/INI); viewer.hpp, polyscope placement (--place)
tests/                     Catch2 v3, one file per module, tagged "[module]"; scikit-fmm fixtures in tests/data/
tools/                     uv single-file Python scripts (PEP 723 + lock file)
assets/meshes/             open bases hemisphere.obj and saddle.obj (radius 1), element ball.obj (radius 0.1)
```

Make `pavel_core` STATIC when the first `.cpp` lands.

## Dependencies

| Dependency | Used for | Paper equivalent |
|---|---|---|
| Eigen | linear algebra, the `V`/`F` mesh matrices | numpy |
| libigl | mesh IO, fast winding number, marching cubes, signed distance | trimesh, MC |
| Geogram | CVT seeding, CVT remeshing | Geogram |
| geometry-central | stripe patterns (Knöppel et al. 2015), smoothest direction fields | Knöppel stripes |
| CLI11, Catch2 | CLI, tests | — |
| polyscope (optional, `PAVEL_VIEWER`) | manual placement viewer | their UI |

- Fast marching is in-house (`src/fast_marching.hpp`, a port of scikit-fmm 2025.6.23's travel-time marcher, quirks included) for per-voxel speed, obstacle masks and per-element grids.
- The only parallelism is concurrent per-element FMM, with `std::thread` and an atomic work index. No OpenMP (Apple clang lacks it).
- Meshes are libigl-style `Eigen::MatrixXd V, Eigen::MatrixXi F`. Conversion to Geogram lives in `src/mesh_convert.hpp` and happens only at module boundaries.
- libigl is declared first: its Eigen (5.0.1) is reused by geometry-central. Don't add a second Eigen. geometry-central builds without SuiteSparse. Geogram builds shared.

## Pipeline (paper → modules)

**Element frame:** +z is the symmetry axis, z = 0 the footprint plane on the base. Elements are placed with +z along the base normal at the seed.

1. **Seeding (§5)**, `seeding.hpp`: manual (seeds file, format in the README; angle 0 = minimal rotation from +z), isotropic (Geogram CVT, N from coverage = N·A_footprint / A_base, default 1.4), offset surface, stripes. Then perturbation: rotation, element choice, scale jitter ≤ ±10%.
2. **Voxelization (§4)**, `voxel_grid.hpp`, `voxelize.hpp`: one dense grid per element on **one global axis-aligned lattice**, so a grid is an integer offset plus extents and Booleans are index arithmetic. The base is voxelized over each element's grid. Paper grid diagonals are roughly 120–200 voxels.
3. **Overlap resolution + volume recovery (§4)**, `deform.hpp`, `fast_marching.hpp`, `recovery_function.hpp`: shared voxels go to the nearest element center, voxels inside the base are removed. Speed `F(x) = Σ_c ω_c · f_{a,b}(r_c(x)/d)`, d = ∛(element volume), defaults a = 0.1, b = 0.3. Phase A computes Tᵢ concurrently, phase B assigns voxels sequentially. **Invariant: an element never exceeds its original voxel count, and reaches it unless its front runs out.**
4. **Output (§6)**, `output.hpp`: marching cubes per element, merge with the base, Geogram CVT remesh (also the smoothing step). Optional printable shells.

## Verified corrections to the printed equations

- **Eq. 1 is misprinted in both versions** (not C¹, doesn't match its own plot). The version matching the plot and the stated properties:
  - `x ≤ 0` → 0
  - `0 < x ≤ a` → −x²/a² + 2x/a
  - `a < x ≤ (a+b)/2` → 1 − 2(x−a)²/(b−a)²
  - `(a+b)/2 < x ≤ b` → 2(x−b)²/(b−a)²
  - `x > b` → 0
- The TOG text writes `|∇T| = F` and `max(D⁻T, D⁺T, 0)`. Eq. 2's 1/F² and standard upwinding imply `|∇T| = 1/F` and `max(D⁻T, −D⁺T, 0)`, as in scikit-fmm's `travel_time`.

## scikit-fmm behaviour matched (see `tools/fmm_reference.py`)

- Binary φ (−1 inside, +1 outside) puts the zero level set halfway between voxels: at speed 1 the first free neighbour gets T = 0.5.
- Speed 0 means unreachable (NaN in the fixtures), not an error. Masked voxels are obstacles.
- Second order gives smaller T than first order (3.3983 vs 3.5 at 4 voxels from a point source).

## Open interpretations

Items marked (decided) were settled with the maintainer; ask before changing them. Ask about the others once they become relevant.

- **Grid "50% larger"**: the element's bounding box scaled by 1.5 about its center, rounded outward to whole voxels.
- **Element center**: the centroid of its voxels before resolution (decided). Ties go to the lower index.
- **Footprint**: the maximal z-section of the element's voxels at the run's voxel size. Its height is the offset d.
- **Isotropic seeds**: each CVT point goes to its closest base point, with interpolated vertex normals.
- **CVT iterations**: Lloyd then Newton, defaulting to `GEO::remesh_smooth`'s 5 and 30. Both are exposed.
- **Contact surface for r_c**: element i's voxels 6-adjacent to neighbor c after overlap removal (decided), with a speed-1 distance field. The base is not a contact (decided): voxels lost to it count in vᵢ but add no speed.
- **ω_c**: the raw count of voxels i lost inside c. A per-element scale of F changes nothing in round-robin phase B.
- **Phase B**: round-robin, one voxel per element per step (arXiv text).
- **Offset surface**: marching cubes of the base's signed distance (fast winding number sign) at d. The CVT runs on it, N comes from the base area (decided), and seeds go to their closest base point.
- **Several elements, automatic seeding**: the mean footprint area and height set N and d.
- **Stripes**: main stripes follow geometry-central's smoothest line field (decided), the second family is turned by α. `--stripes D L α`: main-stripe distance D, seed spacing L along a main stripe, so the second stripes are L·sin α apart. Seeds in singular faces are dropped. The element's +x runs along the main stripe.
- **Stripe connectivity**: estimated (decided). Seeds within 1.5 gaps and ≤ 20° off a family's direction are linked. A BFS aligns the main directions and 2-colours the seeds (flip along main stripes, keep along second ones). That parity drives `alternate`.
- **Rotation modes**: `none` keeps the seeding's orientation, `random` adds a uniform turn, `alternate` adds 180° by parity. Alternation is stripe-only (decided).
- **Merge before remesh**: plain concatenation, no Boolean union (decided).
- **Remesh point count**: `--remesh-points` defaults to the merged vertex count, so the CVT only smooths.
- **Shells**: the axis is the base's z axis through its bounding-box center (decided), cut at x = c.x and y = c.y. Elements join the patch holding their seed and are never cut (decided). Each patch is remeshed without the base at its merged vertex count.
- **Manual placement UI**: a click on the base adds a seed (projected as in isotropic seeding). The panel edits the newest seed, Undo removes it, and Deform previews in the UI callback. Closing the window runs the full pipeline, with no caching (decided).
- **FMM order**: Eq. 2 is first order, but scikit-fmm defaults to second. First order is the default, `--order` switches.
