"""Tests of the MAFFT alignment of junction reads.

A 2 kb random reference with a tandem duplication of 501-1,500 (1-based):
its junction joins 1,500 back to 501. Reads carry 3 extra bases (GAT) at the
junction, and one reference has a 2 bp microhomology (AC) there.
"""
from __future__ import annotations

import random

import pysam
import pytest

from ardu.junction_align import (
    align_junction,
    expected_junction,
    mafft_available,
    microhomology,
    microhomology_split,
)

pytestmark = pytest.mark.skipif(not mafft_available(), reason="mafft not installed")

ENDS = ((501, "L"), (1_500, "R"))


def _fasta(tmp_path, seq):
    path = tmp_path / "ref.fa"
    path.write_text(f">chr1\n{seq}\n")
    pysam.faidx(str(path))
    return pysam.FastaFile(str(path))


def _reference(microhomology_bases=""):
    rng = random.Random(7)
    ref = [rng.choice("ACGT") for _ in range(2_000)]
    # break any chance homology around both junction ends, then add the wanted one
    ref[1_500 - 1 - len(microhomology_bases)] = "A"
    ref[500 - 1 - len(microhomology_bases)] = "C"
    ref[1_500] = "G"
    ref[500] = "T"
    for i, base in enumerate(microhomology_bases):
        ref[1_500 - len(microhomology_bases) + i] = base
        ref[500 - len(microhomology_bases) + i] = base
    return "".join(ref)


def _junction_reads(ref, inserted="GAT", n=8):
    junction = ref[1_500 - 150:1_500] + inserted + ref[500:500 + 150]
    starts = range(150 - 90, 150 - 10, 10)[:n]
    return {f"r{i}": junction[s:s + 100] for i, s in enumerate(starts)}, junction


def test_microhomology(tmp_path):
    assert microhomology(_fasta(tmp_path, _reference("AC")), "chr1", ENDS) == "AC"
    assert microhomology_split(_fasta(tmp_path, _reference("AC")), "chr1", ENDS) == ("AC", "")
    assert microhomology(_fasta(tmp_path, _reference()), "chr1", ENDS) == ""


def test_inserted_bases_against_expected(tmp_path):
    ref = _reference()
    fasta = _fasta(tmp_path, ref)
    expected = expected_junction(fasta, "chr1", ENDS)
    assert expected[0] == ref[1_350:1_500] + ref[500:650]
    reads, _ = _junction_reads(ref)
    stats = align_junction(reads, {}, expected)["stats"]
    assert stats["inserted_bases"] == "GAT"
    assert (stats["aligned_reads"], stats["junction_depth"]) == (8, 8)
    assert stats["junction_agreement"] == 1


def test_reference_free_junction_from_split_reads(tmp_path):
    reads, junction = _junction_reads(_reference(), inserted="")
    # r0 starts 90 bp before the junction, so its 91st base is the first right of it
    result = align_junction(reads, {"r0": (90, 90)}, None)
    assert result["stats"]["junction_depth"] == 8
    assert result["consensus"].upper() == junction[60:230]
    assert result["consensus_junction"] == 90
    marker, consensus = result["text"].splitlines()[:2]
    assert consensus[marker.index("|")] == junction[150]


def test_reference_free_inserted_bases():
    reads, _ = _junction_reads(_reference())
    # r0: 90 bases of the left part, then GAT, then the right part
    stats = align_junction(reads, {"r0": (90, 93), "r1": (80, 83)}, None)["stats"]
    assert stats["inserted_bases"] == "GAT"
