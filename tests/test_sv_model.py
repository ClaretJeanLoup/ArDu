"""Tests of the haplotype model built from joined breakpoint pairs.

Four breakpoints cut the locus into flank 0, segments 1 to 3 and flank 4.
"""
from __future__ import annotations

from ardu.sv_model import build_model

DUP = ("1-4", "duplication", 0, 3)


def _segments(pieces):
    return [p["seg"] for p in pieces]


def test_deletion_in_one_copy_at_random():
    pieces, notes = build_model(4, [DUP, ("2-3", "deletion", 1, 2)], [1, 2.4, 1.3, 2.7, 1])
    assert _segments(pieces) in ([0, 1, 3, 1, 2, 3, 4], [0, 1, 2, 3, 1, 3, 4])
    assert "chosen at random" in notes[0]


def test_deletion_in_both_copies_when_less_covered_than_flanks():
    pieces, notes = build_model(4, [DUP, ("2-3", "deletion", 1, 2)], [1, 2, 0.5, 2, 1])
    assert _segments(pieces) == [0, 1, 3, 1, 3, 4]
    assert "both copies" in notes[0]


def test_internal_duplication():
    pieces, _ = build_model(4, [DUP, ("2-3", "duplication", 1, 2)], [1, 2, 3, 2, 1])
    assert _segments(pieces) in ([0, 1, 2, 2, 3, 1, 2, 3, 4], [0, 1, 2, 3, 1, 2, 2, 3, 4])
    assert [p["junction"] for p in pieces].count("2-3") == 1


def test_crossing_junctions_not_modelled():
    pieces, notes = build_model(4, [("1-3", "duplication", 0, 2), ("2-4", "deletion", 1, 3)])
    assert _segments(pieces) == [0, 1, 2, 1, 2, 3, 4]
    assert "not modelled" in notes[0]


def test_same_seed_same_model():
    events = [DUP, ("2-3", "deletion", 1, 2)]
    assert build_model(4, events, seed=5) == build_model(4, events, seed=5)
