"""Tests of the structure profile: one row per joined deletion or tandem
duplication, drawn to scale with the duplications that contain it.

Four breakpoints cut the locus into flank 0, segments 1 to 3 and flank 4.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ardu.cli.junctions import (_breakpoint_order, _draw_profile, _read_marks, _read_target,
                                _row_levels, _structure_rows)

BKPS = [1000, 2000, 3000, 4000]


def _segments(depths):
    bounds = [500] + BKPS + [4500]
    labels = ["flank", "1-2", "2-3", "3-4", "flank"]
    return [{"bam": "s", "locus": "gene", "segment": lab, "chrom": "2R", "start": a + 1, "end": b,
             "median_depth": d * 20, "relative_depth": d}
            for lab, a, b, d in zip(labels, bounds, bounds[1:], depths)]


def _ev(pair, sv_type, size=None):
    return {"pair": pair, "sv_type": sv_type, "size": size, "support": 10}


def _levels(joined):
    rows = _structure_rows([1, 2, 3, 4], joined)
    return {e[0]["pair"]: _row_levels(5, e, c) for e, c in rows}


def test_breakpoint_order_from_segments():
    assert _breakpoint_order(_segments([1] * 5)) == [1, 2, 3, 4]


def test_duplication_alone():
    assert _levels([_ev("1-4", "duplication")]) == {"1-4": [1, 2, 2, 2, 1]}


def test_deletion_inside_the_duplication_drops_to_one_copy():
    levels = _levels([_ev("1-4", "duplication"), _ev("2-3", "deletion")])
    assert levels == {"1-4": [1, 2, 2, 2, 1], "2-3": [1, 2, 1, 2, 1]}


def test_two_overlapping_deletions_get_their_own_rows():
    # a heterozygote for two internal deletions that overlap without nesting
    levels = _levels([_ev("1-4", "duplication"), _ev("1-3", "deletion"), _ev("2-4", "deletion")])
    assert levels["1-3"] == [1, 1, 1, 2, 1]
    assert levels["2-4"] == [1, 2, 1, 1, 1]


def test_deletion_on_single_copy_sequence_drops_to_zero():
    assert _levels([_ev("2-3", "deletion")]) == {"2-3": [1, 1, 0, 1, 1]}


def test_profile_figure_with_target_and_mark():
    segments = _segments([1, 1.9, 1.4, 1.9, 1])
    joined = [_ev("1-4", "duplication", 3000), _ev("2-3", "deletion", 1000)]
    target = [("2R", 1200, 1300, "gene"), ("2R", 1400, 1600, "gene"), ("2R", 3200, 3300, "other")]
    fig, ax = plt.subplots()
    _draw_profile(ax, "s", "gene", segments, joined, target, [("2R", 1450, "G119S")])
    texts = [t.get_text() for t in ax.texts]
    assert "dup 3.0 kb" in texts and "Δ 1.0 kb" in texts
    assert "gene" in texts and "G119S" in texts
    assert any(t.startswith("2-3 deletion\nin 1-4") for t in texts)
    # only the rows named after the locus are drawn
    assert len(ax.patches) == 2
    plt.close(fig)


def test_read_target_and_marks(tmp_path):
    bed = tmp_path / "loci.bed"
    bed.write_text("2R\t10\t20\tgene\n\n2R\t30\t40\n")
    mut = tmp_path / "mut.tsv"
    mut.write_text("2R\t15\tG119S\n2R\t35\n")
    assert _read_target(str(bed)) == [("2R", 10, 20, "gene"), ("2R", 30, 40, "")]
    assert _read_marks(str(mut)) == [("2R", 15, "G119S"), ("2R", 35, "")]
