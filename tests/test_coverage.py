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


def test_mutation_indels_and_runs(synthetic, tmp_path):
    (tmp_path / "mut.tsv").write_text("chr1\t12000\thomopolymer\nchr1\t3000\tplain\n")
    run_dir = run_coverage(synthetic, tmp_path, "--mutation", str(tmp_path / "mut.tsv"), "--mutation-indels")
    rows = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(os.path.join(run_dir, "test_mutations.tsv"))}
    hp = dict(kv.split("=") for kv in rows["chr1:12000:homopolymer"].split(";"))
    assert hp["ins"] == "8" and hp["del"] == "0"
    assert hp["runs"] == "T8:12,T9:3,T10:5"
    plain = dict(kv.split("=") for kv in rows["chr1:3000:plain"].split(";"))
    assert plain["ins"] == "0" and plain["del"] == "0"


def test_mutation_output_unchanged_without_flag(synthetic, tmp_path):
    (tmp_path / "mut.tsv").write_text("chr1\t3000\tplain\n")
    run_dir = run_coverage(synthetic, tmp_path, "--mutation", str(tmp_path / "mut.tsv"))
    cell = open(os.path.join(run_dir, "test_mutations.tsv")).read().splitlines()[1].split("\t")[1]
    assert [kv.split("=")[0] for kv in cell.split(";")] == ["A", "T", "C", "G", "depth"]


def test_rcn_interval(synthetic, tmp_path):
    rcn = read_rcn(run_coverage(synthetic, tmp_path, "--rcn-ci", "200"))
    low, high = float(rcn["dup"][5]), float(rcn["dup"][6])
    assert low <= 2.0 <= high
    assert high - low < 0.5
    out = subprocess.run([sys.executable, "-m", "ardu", "parse", "-i",
                          glob.glob(os.path.join(tmp_path, "ArDuRun-test-*", "test_coverage.tsv"))[0],
                          "--RCN", "--RCNlow", "--RCNhigh"], check=True, capture_output=True, text=True)
    parsed = glob.glob(os.path.join(tmp_path, "ArDuRun-test-*", "RCN_RCNlow_RCNhigh.tsv"))
    assert parsed and "_RCNlow" in open(parsed[0]).readline()


def test_no_interval_by_default(synthetic, tmp_path):
    rcn = read_rcn(run_coverage(synthetic, tmp_path))
    assert len(rcn["dup"]) == 5


def test_rerun_gets_new_directory(synthetic, tmp_path):
    first = run_coverage(synthetic, tmp_path)
    second = run_coverage(synthetic, tmp_path)
    assert first != second and second.endswith("_2")
    assert os.path.exists(os.path.join(first, "test_coverage.tsv"))


def test_exclude_zero_depth(synthetic, tmp_path):
    # 4 bases at 20x plus 5 uncovered bases at the contig end (the last read ends at 19,995)
    regions = tmp_path / "gap.bed"
    regions.write_text("chr1\t1001\t4000\tnorm\nchr1\t3001\t3004\tgap\nchr1\t19996\t20000\tgap\n")
    cells = {}
    for mode, extra in (("zeros", []), ("nozeros", ["--exclude-zero-depth"])):
        (tmp_path / mode).mkdir()
        cmd = [sys.executable, "-m", "ardu", "coverage", "-b", str(synthetic["bamlist"]),
               "-r", str(regions), "-n", "norm", "-o", "test", *extra]
        subprocess.run(cmd, cwd=tmp_path / mode, check=True, capture_output=True, text=True)
        cells[mode] = read_rcn(glob.glob(os.path.join(tmp_path, mode, "ArDuRun-test-*"))[0])["gap"]
    mean, _, median, covered, rcn = cells["zeros"]
    assert float(median) == 0.0 and float(rcn) == 0.0 and abs(float(mean) - 80 / 9) < 0.01
    mean, _, median, covered_nz, rcn = cells["nozeros"]
    assert float(median) == 20.0 and float(rcn) == 1.0 and float(mean) == 20.0
    assert covered == covered_nz == "4"
