"""Tests of breakpoint detection with a penalty instead of a set number.

Synthetic normalised depth over 300 kb (Poisson reads around 10x, divided by
10), with steps at 60,000 (x2), 120,000 (back to x1.5) and 200,000 (back to x1).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from ardu.breakpoints import refine_breakpoint, run_ruptures_univariate

STEPS = (60_000, 120_000, 200_000)


def _depth(bin_size=None):
    rng = np.random.default_rng(3)
    pos = np.arange(300_000)
    level = np.select([pos < STEPS[0], pos < STEPS[1], pos < STEPS[2]], [1.0, 2.0, 1.5], 1.0)
    d = pd.DataFrame({"pos": pos, "norm": rng.poisson(10 * level) / 10})
    d["base_pos"] = d["pos"]
    if bin_size:
        d["genomic_bins"] = pd.cut(d["pos"], bins=np.arange(0, 300_000 + bin_size, bin_size),
                                   include_lowest=True)
        d["pos"] = d.groupby("genomic_bins", observed=True)["pos"].transform("mean")
    return d


def test_penalty_finds_the_steps_without_a_number():
    for bin_size in (None, 500):
        d = _depth(bin_size)
        for algo in ("Pelt", "BottomUp"):
            bkps = run_ruptures_univariate(d, algo_name=algo, pen=10)
            found = [int(d.loc[b["row_index"], "base_pos"]) for b in bkps]
            assert len(found) == 3
            assert all(abs(f - s) <= 100 for f, s in zip(found, STEPS))


def test_higher_penalty_keeps_only_the_largest_steps():
    d = _depth(500)
    assert len(run_ruptures_univariate(d, algo_name="Pelt", pen=10)) == 3
    assert len(run_ruptures_univariate(d, algo_name="Pelt", pen=1e6)) == 0


def test_refine_breakpoint():
    x = np.r_[np.ones(400), 2 * np.ones(600)]
    assert refine_breakpoint(x, 0, 1000) == 400
    assert refine_breakpoint(x, 200, 900) == 400
