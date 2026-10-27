"""Synthetic BAM fixtures with known answers for the ArDu tests.

Reference ``chr1`` (20 kb, random sequence) with:
    - uniform 20x coverage of 100 bp reads,
    - a duplicated region 5,001-7,000 at 40x (extra reads),
    - a region 9,001-10,000 where 25% of the reads have MAPQ 0,
    - a tandem duplication junction joining 7,000 back to 5,001: 10 mapped
      reads soft-clipped at each side, plus 6 unmapped junction reads, half
      stored as reverse complement,
    - a homopolymer of 8 T at 12,001-12,008 with reads carrying T8, T9 or T10.
"""
from __future__ import annotations

import random

import pysam
import pytest

LEN = 20_000
READ = 100
HP_START, HP_LEN = 12_001, 8          # 1-based start of the T run


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def _reference() -> str:
    rng = random.Random(1)
    seq = [rng.choice("ACG") for _ in range(LEN)]          # no T, so the run below is unambiguous
    seq[HP_START - 1:HP_START - 1 + HP_LEN] = "T" * HP_LEN
    return "".join(seq)


def _segment(header, name, seq, pos=None, cigar=None, mapq=60, flag=0):
    a = pysam.AlignedSegment(header)
    a.query_name = name
    a.query_sequence = seq
    a.flag = flag
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    if pos is None:
        a.reference_id = -1
        a.reference_start = -1
        a.mapping_quality = 0
    else:
        a.reference_id = 0
        a.reference_start = pos
        a.mapping_quality = mapq
        a.cigarstring = cigar or f"{len(seq)}M"
    return a


@pytest.fixture(scope="session")
def synthetic(tmp_path_factory):
    d = tmp_path_factory.mktemp("synthetic")
    ref = _reference()
    fasta = d / "ref.fa"
    fasta.write_text(">chr1\n" + "\n".join(ref[i:i + 60] for i in range(0, LEN, 60)) + "\n")
    pysam.faidx(str(fasta))

    header = pysam.AlignmentHeader.from_dict({"HD": {"VN": "1.6", "SO": "coordinate"},
                                              "SQ": [{"SN": "chr1", "LN": LEN}]})
    reads, n = [], 0

    def add(seq, pos=None, cigar=None, mapq=60, flag=0):
        nonlocal n
        n += 1
        reads.append(_segment(header, f"r{n}", seq, pos, cigar, mapq, flag))

    # 20x everywhere: a read start every 5 bp
    for p in range(0, LEN - READ, 5):
        if HP_START - 1 - READ <= p <= HP_START - 1 + HP_LEN:     # reads over the run are added below
            continue
        mapq = 0 if 9_000 <= p < 10_000 and (p // 5) % 4 == 0 else 60
        add(ref[p:p + READ], p, mapq=mapq)
    # duplicated region 5,001-7,000 gets 20x more
    for p in range(5_000, 7_000 - READ, 5):
        add(ref[p:p + READ], p)
    # tandem duplication junction: ...6,951-7,000 | 5,001-5,050...
    junction = ref[6_950:7_000] + ref[5_000:5_050]
    for i in range(10):
        add(junction, 6_950, cigar="50M50S")                       # right clip at 7,000
        add(junction, 5_000, cigar="50S50M")                      # left clip at 5,001
    for i in range(6):
        add(junction if i % 2 == 0 else _revcomp(junction), None, flag=4)
    # homopolymer reads, 50 bp anchor on both sides, run lengths 8, 9 and 10
    for length, count in ((8, 12), (9, 3), (10, 5)):
        left = ref[HP_START - 1 - 50:HP_START - 1]
        right = ref[HP_START - 1 + HP_LEN:HP_START - 1 + HP_LEN + 50 - (length - HP_LEN)]
        seq = left + "T" * length + right
        cigar = f"50M{length - HP_LEN}I{len(right) + HP_LEN}M" if length > HP_LEN else f"{len(seq)}M"
        for _ in range(count):
            add(seq, HP_START - 1 - 50, cigar=cigar)

    reads.sort(key=lambda a: (a.reference_id < 0, a.reference_start))
    bam = d / "sample.bam"
    with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
        for a in reads:
            out.write(a)
    pysam.index(str(bam))

    regions = d / "regions.bed"
    regions.write_text("chr1\t1001\t4000\tnorm\nchr1\t5201\t6800\tdup\nchr1\t9001\t10000\tlowmapq\n")
    bamlist = d / "bams.txt"
    bamlist.write_text(f"{bam}\n")
    return dict(dir=d, fasta=fasta, bam=bam, regions=regions, bamlist=bamlist, ref=ref,
                junction=junction)
