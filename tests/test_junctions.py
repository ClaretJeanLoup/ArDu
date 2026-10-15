"""Tests of the junction-read search."""
from __future__ import annotations

from ardu.junctionlib import (
    build_softclip_automaton,
    find_reads_with_multiple_softclips_aho,
    write_junction_reads_fastq,
)


def automaton_for(synthetic):
    ref = synthetic["ref"]
    # right clip of reads ending at 7,000 and left clip of reads starting at 5,001
    clips = {1: [ref[5_000:5_050]], 2: [ref[6_950:7_000]]}
    return build_softclip_automaton(clips)


def test_unmapped_reads_found_in_both_orientations(synthetic):
    hits = find_reads_with_multiple_softclips_aho(
        str(synthetic["bam"]), automaton_for(synthetic), ("chr1", 4_900, 7_100))
    names = {h["read_name"] for h in hits}
    unmapped = [h for h in hits if h["strand"] == "-"]
    assert len(names) == 26          # 20 mapped junction reads and all 6 unmapped ones
    assert len(unmapped) == 3        # the 3 stored as reverse complement need the second strand


def test_unplaced_reads_scanned_only_on_request(synthetic):
    hits = find_reads_with_multiple_softclips_aho(
        str(synthetic["bam"]), automaton_for(synthetic), ("chr1", 4_900, 7_100), include_unplaced=False)
    assert len({h["read_name"] for h in hits}) == 20


def test_fastq_written_without_rereading_bam(synthetic, tmp_path):
    hits = find_reads_with_multiple_softclips_aho(
        str(synthetic["bam"]), automaton_for(synthetic), ("chr1", 4_900, 7_100))
    write_junction_reads_fastq("not/a/bam", hits, str(tmp_path), "sample")
    fq = list(tmp_path.glob("*.fastq"))
    assert len(fq) == 1
    lines = fq[0].read_text().splitlines()
    assert len(lines) == 4 * len(hits)
    assert all(len(lines[i + 1]) == len(lines[i + 3]) for i in range(0, len(lines), 4))
