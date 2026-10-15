"""Tests of ``ardu coverage``: RCN precision, read filters, output directory, mutation counts."""
from __future__ import annotations

import glob
import os
import subprocess
import sys


def run_coverage(synthetic, cwd, *extra):
    cmd = [sys.executable, "-m", "ardu", "coverage", "-b", str(synthetic["bamlist"]),
           "-r", str(synthetic["regions"]), "-n", "norm", "-o", "test", *extra]
    subprocess.run(cmd, cwd=cwd, check=True, capture_output=True, text=True)
    runs = sorted(glob.glob(os.path.join(cwd, "ArDuRun-test-*")))
    return runs[-1]


def read_rcn(run_dir):
    rcn = {}
    for line in open(os.path.join(run_dir, "test_coverage.tsv")):
        if line.startswith("#") or line.startswith("locus"):
            continue
        locus, cell = line.rstrip("\n").split("\t")[:2]
        rcn[locus] = cell.split(";")
    return rcn


def test_rcn_has_four_decimals(synthetic, tmp_path):
    rcn = read_rcn(run_coverage(synthetic, tmp_path))
    assert float(rcn["dup"][4]) == 2.0
    assert float(rcn["lowmapq"][4]) == 0.75
    assert len(rcn["lowmapq"][4].split(".")[1]) >= 2      # 0.75 would print as 0.8 with one decimal


def test_min_mapq_default_and_option(synthetic, tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    default = read_rcn(run_coverage(synthetic, tmp_path / "a"))
    lowered = read_rcn(run_coverage(synthetic, tmp_path / "b", "--min-mapq", "0"))
    assert float(default["lowmapq"][4]) == 0.75            # MAPQ 0 reads dropped by default (MAPQ >= 20)
    assert float(lowered["lowmapq"][4]) == 1.0             # and kept with --min-mapq 0
    assert float(lowered["dup"][4]) == 2.0
