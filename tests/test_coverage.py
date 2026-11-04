"""Tests of ``ardu coverage``: RCN precision, read filters, output directory, mutation counts, coverage QC."""
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
    assert float(rcn["dup"][5]) == 2.0
    assert float(rcn["lowmapq"][5]) == 0.75
    assert len(rcn["lowmapq"][5].split(".")[1]) >= 2      # 0.75 would print as 0.8 with one decimal


def test_min_mapq_default_and_option(synthetic, tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    default = read_rcn(run_coverage(synthetic, tmp_path / "a"))
    lowered = read_rcn(run_coverage(synthetic, tmp_path / "b", "--min-mapq", "0"))
    assert float(default["lowmapq"][5]) == 0.75            # MAPQ 0 reads dropped by default (MAPQ >= 20)
    assert float(lowered["lowmapq"][5]) == 1.0             # and kept with --min-mapq 0
    assert float(lowered["dup"][5]) == 2.0


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
    low, high = float(rcn["dup"][6]), float(rcn["dup"][7])
    assert low <= 2.0 <= high
    assert high - low < 0.5
    out = subprocess.run([sys.executable, "-m", "ardu", "parse", "-i",
                          glob.glob(os.path.join(tmp_path, "ArDuRun-test-*", "test_coverage.tsv"))[0],
                          "--RCN", "--RCNlow", "--RCNhigh"], check=True, capture_output=True, text=True)
    parsed = glob.glob(os.path.join(tmp_path, "ArDuRun-test-*", "RCN_RCNlow_RCNhigh.tsv"))
    assert parsed and "_RCNlow" in open(parsed[0]).readline()


def test_no_interval_by_default(synthetic, tmp_path):
    rcn = read_rcn(run_coverage(synthetic, tmp_path))
    assert len(rcn["dup"]) == 6


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
    mean, _, median, covered, _, rcn = cells["zeros"]
    assert float(median) == 0.0 and float(rcn) == 0.0 and abs(float(mean) - 80 / 9) < 0.01
    mean, _, median, covered_nz, _, rcn = cells["nozeros"]
    assert float(median) == 20.0 and float(rcn) == 1.0 and float(mean) == 20.0
    assert covered == covered_nz == "4"


def _run_gap_locus(synthetic, tmp_path, *extra):
    # locus of two regions: 3,001-3,004 at 20x and the uncovered contig end 19,996-20,000
    regions = tmp_path / "gap.bed"
    regions.write_text("chr1\t1001\t4000\tnorm\nchr1\t3001\t3004\tgap\nchr1\t19996\t20000\tgap\n")
    cmd = [sys.executable, "-m", "ardu", "coverage", "-b", str(synthetic["bamlist"]),
           "-r", str(regions), "-n", "norm", "-o", "test", *extra]
    res = subprocess.run(cmd, cwd=tmp_path, check=True, capture_output=True, text=True)
    return glob.glob(os.path.join(tmp_path, "ArDuRun-test-*"))[0], res.stdout


def test_pct_min_cov_and_region_qc(synthetic, tmp_path):
    run_dir, stdout = _run_gap_locus(synthetic, tmp_path)
    cells = read_rcn(run_dir)
    assert float(cells["gap"][4]) == round(100 * 4 / 9, 1)
    assert float(cells["norm"][4]) == 100.0
    qc = [l.rstrip("\n").split("\t") for l in open(os.path.join(run_dir, "test_region_qc.tsv"))
          if not l.startswith("#")]
    assert qc[0] == ["bam_file", "locus", "region", "length", "pct_zero", "pct_ge3x", "flag"]
    gap = {r[2]: r for r in qc[1:] if r[1] == "gap"}
    assert gap["chr1:19996-20000"][4:] == ["100.0", "0.0", "uncovered"]
    assert gap["chr1:3001-3004"][4:] == ["0.0", "100.0", "."]
    assert "gap: 1 of 2 regions" in stdout
    assert {r[1] for r in qc[1:]} == {"norm", "gap"}      # the normalisation locus is reported too


def test_min_cov_and_warn_zero_options(synthetic, tmp_path):
    # at 21x or more nothing of the gap locus is covered; no region reaches 101% at zero
    run_dir, stdout = _run_gap_locus(synthetic, tmp_path, "--min-cov", "21", "--warn-zero", "101")
    assert float(read_rcn(run_dir)["gap"][4]) == 0.0
    assert "pct_ge21x" in open(os.path.join(run_dir, "test_region_qc.tsv")).read()
    assert "Warning: regions" not in stdout


def test_run_time_in_header(synthetic, tmp_path):
    header = open(os.path.join(run_coverage(synthetic, tmp_path), "test_coverage.tsv")).read().splitlines()
    i = next(k for k, l in enumerate(header) if l.startswith("# Generated on:"))
    assert header[i + 1].startswith("# Run time: ") and header[i + 1].endswith(" seconds")


def test_parse_pct_min_cov_with_and_without_interval(synthetic, tmp_path):
    for extra in ([], ["--rcn-ci", "50"]):
        d = tmp_path / ("ci" if extra else "noci")
        d.mkdir()
        cov = os.path.join(run_coverage(synthetic, d, *extra), "test_coverage.tsv")
        subprocess.run([sys.executable, "-m", "ardu", "parse", "-i", cov, "--PctMinCov", "--RCN"],
                       check=True, capture_output=True, text=True)
        out = glob.glob(os.path.join(os.path.dirname(cov), "PctMinCov_RCN.tsv"))[0]
        rows = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1:] for l in open(out)}
        assert float(rows["dup"][1]) == 2.0 and float(rows["norm"][0]) == 100.0


def test_parse_old_file_without_format_line(tmp_path):
    old = tmp_path / "old_coverage.tsv"
    old.write_text("# ArDu coverage outfile\nlocus\ts1\nace-1\t29.35;12.09;31.0;9380;3.4\n")
    subprocess.run([sys.executable, "-m", "ardu", "parse", "-i", str(old), "--RCN"],
                   check=True, capture_output=True, text=True)
    assert open(tmp_path / "RCN.tsv").read().splitlines()[1].split("\t")[1] == "3.4"
