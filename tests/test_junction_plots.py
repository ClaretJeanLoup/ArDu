"""Tests of the alignment panel of the junction evidence plots."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ardu.cli.junctions import _draw_alignment_panel
from ardu.junction_align import EXPECTED, consensus_columns, junction_stats


def _alignment(aligned: dict[str, str], junction_col: int) -> dict:
    columns = consensus_columns(aligned)
    return {"aligned": aligned, "columns": columns, "junction_col": junction_col,
            "stats": junction_stats(aligned, columns, junction_col - 1, junction_col)}


def test_reads_away_from_the_junction():
    # Adapter dimers aligned at the start, the expected junction 100 columns on.
    expected = "ACGT" * 40
    read = "GATCGGAAGAGC" + "-" * (len(expected) - 12)
    alignment = _alignment({EXPECTED: expected, "r1": read, "r2": read}, junction_col=100)
    fig, ax = plt.subplots()
    _draw_alignment_panel(ax, alignment)
    assert ax.images[0].get_array().shape == (3, 12)
    assert not ax.lines
    plt.close(fig)


def test_junction_marked_when_on_the_reads():
    expected = "ACGT" * 40
    read = "-" * 80 + expected[80:]
    alignment = _alignment({EXPECTED: expected, "r1": read}, junction_col=100)
    fig, ax = plt.subplots()
    _draw_alignment_panel(ax, alignment)
    assert ax.images[0].get_array().shape == (2, 50)
    assert len(ax.lines) == 1
    plt.close(fig)
