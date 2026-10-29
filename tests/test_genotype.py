"""Tests of ``ardu genotype``: codon reading, copy sharing and the output table.

Reference ``chr1`` (2 kb, random A/C/G) with two codons:
    - ``fwd`` at 1,001-1,003, GGC (G) on the forward strand, read as
      20 GGC + 10 AGC (S) + 1 TGC (C, a lone error),
    - ``rev`` at 1,501-1,503, GCC on the forward strand, so GGC (G) on the
      reverse one, read as 12 GCC + 12 GCT (AGC, S).
Plus, over ``fwd``: an overlapping pair whose two mates read AGC (counted
once), a read with low quality at the codon, a MAPQ 10 read and a read with a
deletion in the codon, all to be left out.
"""
from __future__ import annotations

import glob
import os
import random
import subprocess
import sys

import pysam
import pytest

from ardu.codons import (
    allocate_copies,
    genotype_string,
    parse_mutation_name,
    total_copies,
)

LEN, READ = 2_000, 100
FWD, REV = 1_001, 1_501               # 1-based first base of each codon


def _reference() -> str:
    rng = random.Random(2)
    seq = [rng.choice("ACG") for _ in range(LEN)]
    seq[FWD - 1:FWD + 2] = "GGC"
    seq[REV - 1:REV + 2] = "GCC"
    return "".join(seq)


def _read(header, name, seq, start, flag=0, mapq=60, quals=None, cigar=None):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.query_sequence = seq
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar or f"{len(seq)}M"
    a.query_qualities = pysam.qualitystring_to_array(quals or "I" * len(seq))
    return a


@pytest.fixture(scope="module")
def codon_data(tmp_path_factory):
    d = tmp_path_factory.mktemp("codons")
    ref = _reference()
    fasta = d / "ref.fa"
    fasta.write_text(">chr1\n" + ref + "\n")
    pysam.faidx(str(fasta))
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "chr1", "LN": LEN}]})
    reads, n = [], 0

    def with_codon(pos, codon):
        start = pos - 51
        seq = ref[start:pos - 1] + codon + ref[pos + 2:start + READ]
        return start, seq

    def add(seq, start, **kw):
        nonlocal n
        n += 1
        reads.append(_read(header, f"r{n}", seq, start, **kw))

    for pos, codon, count in ((FWD, "GGC", 20), (FWD, "AGC", 10), (FWD, "TGC", 1),
                              (REV, "GCC", 12), (REV, "GCT", 12)):
        start, seq = with_codon(pos, codon)
        for _ in range(count):
            add(seq, start)

    # an overlapping pair, both mates AGC: one codon
    start, seq = with_codon(FWD, "AGC")
    m1 = _read(header, "pair", seq, start, flag=1 | 2 | 32 | 64)
    m2 = _read(header, "pair", ref[start + 30:FWD - 1] + "AGC" + ref[FWD + 2:start + 30 + READ],
               start + 30, flag=1 | 2 | 16 | 128)
    for a, mate in ((m1, m2), (m2, m1)):
        a.next_reference_id, a.next_reference_start = 0, mate.reference_start
    reads += [m1, m2]
    # left out: base quality 10 at the codon, MAPQ 10, a deletion of the codon's middle base
    start, seq = with_codon(FWD, "AGC")
    add(seq, start, quals="I" * 50 + "+++" + "I" * 47)
    add(seq, start, mapq=10)
    add(seq[:51] + seq[52:], start, cigar="51M1D48M")

    reads.sort(key=lambda a: a.reference_start)
    bam = d / "sample.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(str(bam))
    (d / "bams.txt").write_text(f"{bam}\n")
    (d / "mutations.tsv").write_text(
        f"chr1\t{FWD}\t{FWD + 2}\tgeneA\tG119S\tforward\n"
        f"chr1\t{REV}\t{REV + 2}\tgeneB\tG10S/A\t-\n"
    )
    (d / "coverage.tsv").write_text(
        "# ArDu coverage outfile\nlocus\tsample\n"
        "geneA\t60;5;60;100;2.0\ngeneB\t30;5;30;100;1.2\n"
    )
    return d


def run_genotype(d, cwd, *extra, check=True):
    cmd = [sys.executable, "-m", "ardu", "genotype", "-b", str(d / "bams.txt"),
           "-m", str(d / "mutations.tsv"), "-f", str(d / "ref.fa"), "-o", "test", *extra]
    return subprocess.run(cmd, cwd=cwd, check=check, capture_output=True, text=True)


def read_cells(cwd):
    cells = {}
    for line in open(os.path.join(cwd, "test_codons.tsv")):
        if line.startswith("#") or line.startswith("mutation"):
            continue
        mid, cell = line.rstrip("\n").split("\t")
        cells[mid.split(":")[0]] = cell.split(";")
    return cells


def test_parse_mutation_name():
    assert parse_mutation_name("G119S") == ("G", frozenset("S"))
    assert parse_mutation_name("L1014S/F") == ("L", frozenset("SF"))
    with pytest.raises(ValueError):
        parse_mutation_name("kdr")


def test_total_copies():
    assert total_copies(None, 2, 2, 1.4) == 2
    assert total_copies(1.3, 2, 2, 1.4) == 2          # noise, not a duplication
    assert total_copies(2.0, 2, 2, 1.4) == 4
    assert total_copies(0.5, 1, 2, 1.4) == 1          # X in a male, normalised on an autosome
    assert total_copies(1.0, 1, 2, 1.4) == 2          # duplicated X in a male
    assert total_copies(3.0, 1, 1, 1.4) == 3          # haploid genome


def test_allocate_copies_sums_and_threshold():
    copies = allocate_copies({"G": 50, "S": 30, "A": 20}, 3, min_reads=1)
    assert copies == {"G": 1, "S": 1, "A": 1}
    assert sum(allocate_copies({"G": 61, "S": 39}, 5).values()) == 5
    # 5% of the reads: under half a copy with 4 copies (12.5%), a copy with 10 (5%)
    assert "S" not in allocate_copies({"G": 95, "S": 5}, 4)
    assert allocate_copies({"G": 94, "S": 6}, 10) == {"G": 9, "S": 1}
    # at least min_reads reads
    assert "S" not in allocate_copies({"G": 2, "S": 2}, 2, min_reads=3)


def test_genotype_string():
    assert genotype_string({"S": 2, "R": 1, "O": 0}) == "2S1R"
    assert genotype_string({"S": 0, "R": 1, "O": 0}) == "1R"


def test_codons_and_genotypes_without_coverage(codon_data, tmp_path):
    run_genotype(codon_data, tmp_path)
    cells = read_cells(tmp_path)
    genotype, s, r, o, depth, rcn, codons, flags = cells["G119S"]
    # 20 GGC, 10 AGC + the pair once; low quality, MAPQ 10 and deletion reads left out
    assert (s, r, o, depth) == ("20", "11", "1", "32")
    assert codons == "GGC:20:G,AGC:11:S,TGC:1:C"
    assert genotype == "1S1R" and rcn == "NA"
    assert flags == "indel_reads_1"
    genotype, s, r, o, depth, *_ = cells["G10S/A"]
    assert (s, r, depth) == ("12", "12", "24") and genotype == "1S1R"


def test_genotypes_with_coverage(codon_data, tmp_path):
    run_genotype(codon_data, tmp_path, "-c", str(codon_data / "coverage.tsv"))
    cells = read_cells(tmp_path)
    assert cells["G119S"][0] == "3S1R" and cells["G119S"][5] == "2"     # 4 copies, 20:11
    assert cells["G10S/A"][0] == "1S1R"                                 # RCN 1.2, not duplicated


def test_haploid_and_low_depth(codon_data, tmp_path):
    run_genotype(codon_data, tmp_path, "--haploid", "all", "--min-depth", "30")
    cells = read_cells(tmp_path)
    assert cells["G119S"][0] == "1S"
    assert cells["G10S/A"][0] == "NA" and "low_depth" in cells["G10S/A"][7]


def test_reference_mismatch_stops(codon_data, tmp_path):
    bad = tmp_path / "bad.tsv"
    bad.write_text(f"chr1\t{FWD}\t{FWD + 2}\tgeneA\tA119S\tforward\n")
    res = subprocess.run([sys.executable, "-m", "ardu", "genotype", "-b", str(codon_data / "bams.txt"),
                          "-m", str(bad), "-f", str(codon_data / "ref.fa"), "-o", "test"],
                         cwd=tmp_path, capture_output=True, text=True)
    assert res.returncode != 0 and "gives G, not A" in res.stderr
    assert not glob.glob(str(tmp_path / "test_codons.tsv"))
