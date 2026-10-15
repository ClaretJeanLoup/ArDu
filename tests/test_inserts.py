"""Tests of the insert size and orientation of read pairs in ``ardu junctions``.

Reference ``chr1`` (10 kb) with 100 bp reads in FR pairs, fragments of
280-320 bp (median 300, MAD 10), and two loci both numbered 1 and 2:
    - ``del``: a 600 bp deletion between 2,000 and 2,600, spanned by 10 FR pairs
      900 bp apart,
    - ``dup``: a 2 kb tandem duplication between 5,000 and 7,000, whose junction
      is crossed by 10 RF pairs,
plus a duplicate and a MAPQ 0 copy of a deletion pair, which must be ignored,
and 5 split reads at each junction: 2,000 | 2,601 for the deletion and
7,000 | 5,001 for the duplication (1-based).
"""
from __future__ import annotations

import argparse

import pandas as pd
import pysam
import pytest

from ardu.bam import read_pair_info
from ardu.cli import junctions
from ardu.cli.junctions import _clipped_end, _sv_geometry

LEN, READ = 10_000, 100
PAIRED, PROPER, REVERSE, MATE_REVERSE, READ1, READ2, DUP = 1, 2, 16, 32, 64, 128, 1024


def _pair(header, name, fwd_start, rev_start, proper=True, extra_flag=0, mapq=60):
    """A forward mate at *fwd_start* and a reverse mate at *rev_start* (0-based)."""
    tlen = abs(rev_start + READ - fwd_start)
    mates = []
    for start, mate_start, flag in (
        (fwd_start, rev_start, PAIRED | MATE_REVERSE | READ1),
        (rev_start, fwd_start, PAIRED | REVERSE | READ2),
    ):
        a = pysam.AlignedSegment(header)
        a.query_name = name
        a.query_sequence = "A" * READ
        a.query_qualities = pysam.qualitystring_to_array("I" * READ)
        a.flag = flag | (PROPER if proper else 0) | extra_flag
        a.reference_id = a.next_reference_id = 0
        a.reference_start, a.next_reference_start = start, mate_start
        a.mapping_quality = mapq
        a.cigarstring = f"{READ}M"
        a.template_length = tlen if start <= mate_start else -tlen
        a.set_tag("MC", f"{READ}M")
        mates.append(a)
    return mates


def _split(header, name, start, cigar, sa_pos, sa_cigar):
    """An unpaired read at 0-based *start* whose clipped part aligns at 1-based *sa_pos*."""
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.query_sequence = "C" * READ
    a.query_qualities = pysam.qualitystring_to_array("I" * READ)
    a.reference_id, a.reference_start = 0, start
    a.mapping_quality = 60
    a.cigarstring = cigar
    a.set_tag("SA", f"chr1,{sa_pos},+,{sa_cigar},60,0;")
    return a


@pytest.fixture(scope="module")
def paired(tmp_path_factory):
    d = tmp_path_factory.mktemp("paired")
    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "chr1", "LN": LEN}]})
    reads = []
    for i, p in enumerate(range(500, LEN - 600, 10)):
        fragment = 280 + (i % 5) * 10
        reads += _pair(header, f"c{i}", p, p + fragment - READ)
    for i in range(10):
        a = 20 + 15 * i                       # bases of the fragment left of the breakpoint
        reads += _pair(header, f"del{i}", 2_000 - a, 2_600 + 300 - a - READ, proper=False)
        reads += _pair(header, f"dup{i}", 7_000 - a, 5_000 + 300 - a - READ, proper=False)
    reads += _pair(header, "del_dup", 1_950, 2_750, proper=False, extra_flag=DUP)
    reads += _pair(header, "del_mapq0", 1_950, 2_750, proper=False, mapq=0)
    for i in range(5):
        reads.append(_split(header, f"sdel{i}", 1_940, "60M40S", 2_601, "60H40M"))
        reads.append(_split(header, f"sdup{i}", 6_950, "50M50S", 5_001, "50H50M"))

    reads.sort(key=lambda a: a.reference_start)
    bam = d / "sample.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(str(bam))

    (d / "bams.txt").write_text(f"{bam}\n")
    (d / "bkp.tsv").write_text(
        "bam_file\tlocus\tchromosome\tposition\tbreakpoint_number\tmethod\n"
        "sample\tdel\tchr1\t2000\t1\truptures\n"
        "sample\tdel\tchr1\t2600\t2\truptures\n"
        "sample\tdup\tchr1\t5000\t1\truptures\n"
        "sample\tdup\tchr1\t7000\t2\truptures\n"
    )
    parser = argparse.ArgumentParser()
    junctions.register(parser.add_subparsers())
    args = parser.parse_args([
        "junctions", "-b", str(d / "bams.txt"), "-i", str(d / "bkp.tsv"),
        "-o", str(d / "out"), "-e", "200", "--pairs", "1-2", "--plot-insert",
        "--plot-softclip", "--plot-evidence",
    ])
    args.func(args)
    return dict(dir=d, bam=bam, header=header,
                summary=pd.read_csv(d / "out_insert_summary.tsv", sep="\t"),
                bridges=pd.read_csv(d / "out_bridging_pairs.tsv", sep="\t"),
                evidence=pd.read_csv(d / "out_junction_summary.tsv", sep="\t").set_index("locus"))


def test_pair_info_fr_and_rf(paired):
    fwd, rev = _pair(paired["header"], "x", 1_000, 1_200)
    for read in (fwd, rev):
        info = read_pair_info(read)
        assert (info.insert, info.orientation) == (300, "FR")
    fwd, rev = _pair(paired["header"], "y", 3_000, 1_000)
    assert read_pair_info(fwd).orientation == read_pair_info(rev).orientation == "RF"


def test_pair_info_filters(paired):
    fwd, _ = _pair(paired["header"], "x", 1_000, 1_200)
    fwd.flag |= 2048                          # supplementary
    assert read_pair_info(fwd) is None
    fwd, _ = _pair(paired["header"], "x", 1_000, 1_200, extra_flag=DUP)
    assert read_pair_info(fwd) is None


def test_library_insert(paired):
    lib = paired["summary"].query("`class` == 'library'").iloc[0]
    assert (lib["median"], lib["mad"]) == (300, 10)


def test_deletion_bridging_pairs(paired):
    row = paired["summary"].query("locus == 'del' and `class` == 'bridging'")
    assert list(row["orientation"]) == ["FR"]
    row = row.iloc[0]
    assert (row["n"], row["median"], row["sv_type"], row["expected"]) == (10, 900, "deletion", 900)


def test_duplication_bridging_pairs(paired):
    row = paired["summary"].query("locus == 'dup' and `class` == 'bridging'")
    assert list(row["orientation"]) == ["RF"]
    row = row.iloc[0]
    assert (row["n"], row["sv_type"], row["expected"]) == (10, "duplication", 1700)
    assert abs(row["median"] - row["expected"]) <= 2


def test_bridging_pairs_listed_once_without_filtered_reads(paired):
    names = paired["bridges"]["read_name"]
    assert names.is_unique
    assert len(names) == 20
    assert not names.isin(["del_dup", "del_mapq0"]).any()
    assert (paired["dir"] / "out_sample_insert_dist.png").exists()


def test_clipped_end_and_geometry():
    assert _clipped_end(1_941, "60M40S") == (2_000, "R")
    assert _clipped_end(2_601, "60H40M") == (2_601, "L")
    assert _clipped_end(100, "100M") is None
    assert _sv_geometry((2_000, "R"), (2_601, "L")) == ("deletion", 600)
    assert _sv_geometry((7_000, "R"), (5_001, "L")) == ("duplication", 2_000)


def test_junction_summary(paired):
    ev = paired["evidence"]
    for locus, b1, b2, sv_type, size, expected in (
        ("del", 2_000, 2_601, "deletion", 600, 900),
        ("dup", 5_001, 7_000, "duplication", 2_000, 1_700),
    ):
        row = ev.loc[locus]
        assert (row["bkp1_refined"], row["bkp2_refined"]) == (b1, b2)
        assert (row["refined_by"], row["sv_type"], row["size"]) == ("split_reads", sv_type, size)
        assert (row["split_reads"], row["bridging_pairs"]) == (5, 10)
        assert row["insert_expected"] == expected
        assert row["junction_reads"] == "NA" or pd.isna(row["junction_reads"])
    for plot in ("softclip_dist", "del_1-2_evidence", "dup_1-2_evidence"):
        assert (paired["dir"] / f"out_sample_{plot}.png").exists()
