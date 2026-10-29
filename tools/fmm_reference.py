# /// script
# requires-python = ">=3.12"
# dependencies = ["numpy==2.5.4", "scikit-fmm==2025.6.23"]
# ///
"""Reference arrival times from scikit-fmm, the FMM the PAVEL authors used (§4).

Writes one directory per case to tests/data/fmm/ holding the inputs (phi, speed, mask)
and the arrival times of first- and second-order FMM (t1, t2). All arrays are C-order
.npy, indexed [i, j, k] with k fastest, dx = 1; NaN marks unreachable voxels.

Run from the repo root: uv run tools/fmm_reference.py
"""

from pathlib import Path

import numpy as np
import skfmm

N = 32
C = N // 2
OUT = Path(__file__).resolve().parents[1] / "tests" / "data" / "fmm"


def point_source():
    phi = np.ones((N, N, N))
    phi[C, C, C] = -1
    return phi, np.ones_like(phi), np.zeros(phi.shape, bool)


def ball_variable_speed():
    i, j, k = np.indices((N, N, N), dtype=float) - C
    phi = np.where(i**2 + j**2 + k**2 <= 36, -1.0, 1.0)
    speed = 1 + 0.5 * np.sin(i / 3) * np.cos(j / 5)
    return phi, speed, np.zeros(phi.shape, bool)


def wall_and_dead_zone():
    phi, speed, mask = point_source()
    mask[:, :, C + 4] = True
    mask[C - 2 : C + 3, C - 2 : C + 3, C + 4] = False
    speed[: N // 4] = 0
    return phi, speed, mask


for case in (point_source, ball_variable_speed, wall_and_dead_zone):
    phi, speed, mask = case()
    out = OUT / case.__name__
    out.mkdir(parents=True, exist_ok=True)
    for name, array in {"phi": phi, "speed": speed, "mask": mask}.items():
        np.save(out / f"{name}.npy", array)
    for order in (1, 2):
        t = np.ma.filled(skfmm.travel_time(np.ma.MaskedArray(phi, mask), speed, order=order), np.nan)
        np.save(out / f"t{order}.npy", t)
        print(f"{case.__name__} order {order}: reachable {np.isfinite(t).sum()}, max T {np.nanmax(t):.4f}")
