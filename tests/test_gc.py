"""Tests of the GC check and of normalisation on random windows."""
from __future__ import annotations

import glob
import os
import subprocess
import sys

import numpy as np
import pysam
import pytest

from ardu.gc import draw_norm_windows, region_gc, write_gc_report
from test_coverage import read_rcn

RANDOM_ARGS = ("--norm-random", "12", "--norm-size", "500", "--norm-margin", "500",
               "--norm-chroms", "chr1")


def run(synthetic, cwd, *extra, norm=("-n", "norm")):
    cmd = [sys.executable, "-m", "ardu", "coverage", "-b", str(synthetic["bamlist"]),
           "-r", str(synthetic["regions"]), *norm, "-o", "test", *extra]
    out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    runs = sorted(glob.glob(os.path.join(cwd, "ArDuRun-test-*")))
    return out, runs[-1] if runs else None


@pytest.fixture
def gc_fasta(tmp_path):
    """One contig of 30 kb: 10 kb at 30% GC, 10 kb at 50%, 10 kb at 70%."""
    rng = np.random.default_rng(5)
    parts = []
    for gc in (0.3, 0.5, 0.7):
        p = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        parts.append("".join(rng.choice(list("ACGT"), size=10_000, p=p)))
    fa = tmp_path / "gc.fa"
    fa.write_text(">c1\n" + "".join(parts) + "\n")
    pysam.faidx(str(fa))
    return str(fa)


def test_region_gc_pads_and_merges(gc_fasta):
    with pysam.FastaFile(gc_fasta) as fa:
        gc, n = region_gc(fa, [("c1", 1_000, 2_000)], pad=0)
        assert n == 1_000 and abs(gc - 0.3) < 0.04
        # two regions whose padding overlaps are counted once
        _, n = region_gc(fa, [("c1", 1_000, 1_100), ("c1", 1_200, 1_300)], pad=100)
        assert n == 500
        # padding is clipped at the contig ends
        _, n = region_gc(fa, [("c1", 0, 100)], pad=250)
        assert n == 350


def test_windows_follow_the_target_gc(gc_fasta):
    w = draw_norm_windows(gc_fasta, 15, 500, exclude=[], target_gc=0.7, chroms=["c1"], pad=0, tolerance=0.05)
    assert len(w) == 15
    assert all(s >= 20_000 for _, s, _ in w)               # only the 70% third matches
    assert w == sorted(w)
    assert all(e - s == 500 for _, s, e in w)
    assert len({s for _, s, _ in w}) == 15                  # no window drawn twice


def test_windows_avoid_excluded_regions_and_are_reproducible(gc_fasta):
    exclude = [("c1", 0, 25_000)]
    a = draw_norm_windows(gc_fasta, 8, 500, exclude, chroms=["c1"], seed=3)
    b = draw_norm_windows(gc_fasta, 8, 500, exclude, chroms=["c1"], seed=3)
    assert a == b and len(a) == 8
    assert all(s >= 25_000 for _, s, _ in a)


def test_too_few_matching_windows_returns_what_exists(gc_fasta):
    w = draw_norm_windows(gc_fasta, 50, 1_000, exclude=[], target_gc=0.7, chroms=["c1"], pad=0)
    assert 0 < len(w) <= 10


def test_gc_report_flags_large_differences(tmp_path):
    flagged = write_gc_report(str(tmp_path / "gc.tsv"),
                              {"ref": (0.55, 1000), "close": (0.52, 900), "far": (0.45, 800)},
                              "ref", warn_points=5)
    assert flagged == ["far"]
    rows = [l.split("\t") for l in open(tmp_path / "gc.tsv") if not l.startswith("#")]
    assert rows[0][:3] == ["locus", "bases", "gc"]
    assert rows[3][3] == "-10.0" and rows[3][4].strip() == "yes"


def test_random_windows_give_the_same_rcn(synthetic, tmp_path):
    out, run_dir = run(synthetic, tmp_path, "--fasta", str(synthetic["fasta"]), *RANDOM_ARGS[2:],
                       norm=RANDOM_ARGS[:2])
    assert out.returncode == 0, out.stderr
    rcn = read_rcn(run_dir)
    assert abs(float(rcn["dup"][4]) - 2.0) < 0.05
    assert abs(float(rcn["norm"][4]) - 1.0) < 0.05
    bed = [l.split("\t") for l in open(os.path.join(run_dir, "test_norm_windows.bed"))]
    assert len(bed) == 12 and all(f[3].strip() == "norm_random" for f in bed)
    # windows keep --norm-margin (500 bp) away from every locus of the region file
    loci = [(1_000, 4_000), (5_200, 6_800), (9_000, 10_000)]
    assert all(int(f[2]) <= s - 500 or int(f[1]) - 1 >= e + 500 for f in bed for s, e in loci)
    assert "norm_random" in open(os.path.join(run_dir, "test_gc.tsv")).read()


def test_windows_reused_with_n(synthetic, tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    _, first = run(synthetic, tmp_path / "a", "--fasta", str(synthetic["fasta"]), *RANDOM_ARGS[2:],
                   norm=RANDOM_ARGS[:2])
    bed = os.path.join(first, "test_norm_windows.bed")
    out, second = run(synthetic, tmp_path / "b", norm=("-n", bed))
    assert out.returncode == 0, out.stderr
    assert read_rcn(first)["dup"] == read_rcn(second)["dup"]


def test_norm_random_needs_fasta_and_excludes_n(synthetic, tmp_path):
    out, _ = run(synthetic, tmp_path, *RANDOM_ARGS[2:], norm=RANDOM_ARGS[:2])
    assert out.returncode != 0 and "--fasta" in out.stderr
    out, _ = run(synthetic, tmp_path, "--fasta", str(synthetic["fasta"]), *RANDOM_ARGS)
    assert out.returncode != 0 and "not allowed with" in out.stderr


def test_gc_report_with_a_named_locus(synthetic, tmp_path):
    out, run_dir = run(synthetic, tmp_path, "--fasta", str(synthetic["fasta"]))
    assert out.returncode == 0, out.stderr
    rows = {l.split("\t")[0]: l.rstrip("\n").split("\t")
            for l in open(os.path.join(run_dir, "test_gc.tsv")) if not l.startswith("#")}
    assert set(rows) == {"locus", "norm", "dup", "lowmapq"}
    assert rows["norm"][3] == "0.0"
    out, run_dir = run(synthetic, tmp_path)
    assert "GC of the loci not checked" in out.stdout
