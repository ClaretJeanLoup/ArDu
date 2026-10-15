"""
Alignment of the reads crossing an SV junction with MAFFT, their consensus,
and the junction sequence expected from the reference.

Covers:
    - expected junction sequence and microhomology from the reference
    - MAFFT alignment of the junction and split reads (L-INS-i, or
      ``--addfragments`` onto the expected junction when there is one)
    - per-column consensus, junction statistics and a text view
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from collections import Counter

FLANK = 150
EXPECTED = "expected"


def mafft_available() -> bool:
    return shutil.which("mafft") is not None


def _junction_ends(ends) -> tuple[int, int] | None:
    """``(right_clipped_end, left_clipped_end)`` of a junction, 1-based.

    The junction joins the reference up to the right-clipped end with the
    reference from the left-clipped end (a deletion or a tandem duplication).
    """
    sides = {side: pos for pos, side in ends}
    if set(sides) != {"L", "R"}:
        return None
    return sides["R"], sides["L"]


def expected_junction(fasta, chrom: str, ends, flank: int = FLANK) -> tuple[str, int] | None:
    """Reference sequence across the junction and the length of its left part."""
    junction = _junction_ends(ends)
    if junction is None:
        return None
    r_end, l_start = junction
    left = fasta.fetch(chrom, max(0, r_end - flank), r_end).upper()
    right = fasta.fetch(chrom, l_start - 1, l_start - 1 + flank).upper()
    return left + right, len(left)


def microhomology(fasta, chrom: str, ends, max_len: int = 50) -> str | None:
    """Bases found on both sides of the junction, so that the breakpoint could
    sit anywhere within them; ``""`` when there are none.
    """
    junction = _junction_ends(ends)
    if junction is None:
        return None
    r_end, l_start = junction
    m = min(max_len, r_end, l_start - 1)
    left = fasta.fetch(chrom, r_end - m, r_end + m).upper()
    right = fasta.fetch(chrom, l_start - 1 - m, l_start - 1 + m).upper()
    k_left = 0
    while k_left < m and left[m - 1 - k_left] == right[m - 1 - k_left]:
        k_left += 1
    k_right = 0
    while k_right < min(len(left), len(right)) - m and left[m + k_right] == right[m + k_right]:
        k_right += 1
    return right[m - k_left:m] + left[m:m + k_right]


def run_mafft(reads: dict[str, str], expected: str | None = None) -> dict[str, str]:
    """Align *reads* (``{name: sequence}``, same strand) and return the aligned
    rows, with the expected junction as row ``EXPECTED`` when given.
    """
    names = list(reads)
    with tempfile.TemporaryDirectory() as tmp:
        reads_fa = os.path.join(tmp, "reads.fa")
        with open(reads_fa, "w") as fh:
            for i, name in enumerate(names):
                fh.write(f">r{i}\n{reads[name]}\n")
        cmd = ["mafft", "--quiet", "--thread", "1"]
        if expected:
            backbone = os.path.join(tmp, "expected.fa")
            with open(backbone, "w") as fh:
                fh.write(f">{EXPECTED}\n{expected}\n")
            cmd += ["--addfragments", reads_fa, "--multipair", backbone]
        else:
            cmd += ["--localpair", "--maxiterate", "1000", reads_fa]
        out = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout

    aligned: dict[str, str] = {}
    for record in out.split(">")[1:]:
        header, *lines = record.split("\n")
        seq = "".join(lines).upper()
        aligned[EXPECTED if header == EXPECTED else names[int(header[1:])]] = seq
    return aligned


def _span(row: str) -> tuple[int, int]:
    stripped = row.strip("-")
    start = len(row) - len(row.lstrip("-"))
    return start, start + len(stripped) - 1


def _column_of_base(row: str, index: int) -> int | None:
    """Alignment column of the *index*-th (0-based) base of *row*."""
    n = -1
    for col, ch in enumerate(row):
        if ch != "-":
            n += 1
            if n == index:
                return col
    return None


def consensus_columns(aligned: dict[str, str]) -> list[tuple[str, int, float]]:
    """``(symbol, depth, agreement)`` per column over the reads (not the
    expected row). A read counts only between its first and last base, so its
    end gaps are not evidence of a deletion.
    """
    rows = [row for name, row in aligned.items() if name != EXPECTED]
    spans = [_span(row) for row in rows]
    columns = []
    for c in range(len(rows[0])):
        counts = Counter(row[c] for row, (a, b) in zip(rows, spans) if a <= c <= b)
        depth = sum(counts.values())
        if not depth:
            columns.append(("", 0, 0.0))
            continue
        symbol, n = counts.most_common(1)[0]
        columns.append((symbol, depth, n / depth))
    return columns


def consensus_sequence(columns) -> str:
    """Consensus without gaps, upper case where at least 3 reads agree at 80%."""
    return "".join(
        symbol if depth >= 3 and agreement >= 0.8 else symbol.lower()
        for symbol, depth, agreement in columns
        if symbol not in ("", "-")
    )


def junction_column(aligned: dict[str, str], expected_offset: int | None,
                    split_offsets: dict[str, int]) -> int | None:
    """Column of the first base right of the junction: from the expected row
    when there is one, else the column most split reads put it at.
    """
    if EXPECTED in aligned and expected_offset is not None:
        return _column_of_base(aligned[EXPECTED], expected_offset)
    votes = Counter(
        col for name, offset in split_offsets.items()
        if name in aligned and (col := _column_of_base(aligned[name], offset)) is not None
    )
    return votes.most_common(1)[0][0] if votes else None


def junction_stats(aligned: dict[str, str], columns, junction_col: int | None,
                   expected_offset: int | None, window: int = 10) -> dict:
    """Depth and lowest agreement around the junction, inserted bases and
    identity of the consensus with the expected junction.
    """
    stats = {"aligned_reads": sum(name != EXPECTED for name in aligned),
             "consensus_length": len(consensus_sequence(columns)),
             "junction_depth": None, "junction_agreement": None,
             "inserted_bases": None, "expected_identity": None}
    if junction_col is not None:
        stats["junction_depth"] = columns[junction_col][1]
        near = [agreement for symbol, depth, agreement in
                columns[max(0, junction_col - window):junction_col + window] if depth]
        stats["junction_agreement"] = round(min(near), 3) if near else None

    expected = aligned.get(EXPECTED)
    if expected is None:
        return stats
    left_col = _column_of_base(expected, expected_offset - 1)
    if left_col is not None and junction_col is not None:
        stats["inserted_bases"] = "".join(
            columns[c][0] for c in range(left_col + 1, junction_col)
            if expected[c] == "-" and columns[c][0] not in ("", "-")
        )
    compared = matches = 0
    for c, (symbol, depth, _) in enumerate(columns):
        if not depth or (symbol == "-" and expected[c] == "-"):
            continue
        compared += 1
        matches += symbol == expected[c]
    stats["expected_identity"] = round(matches / compared, 4) if compared else None
    return stats


def alignment_text(aligned: dict[str, str], columns, junction_col: int | None) -> str:
    """Reads stacked under the consensus, ``.`` where a read agrees with it,
    ``|`` marking the first base right of the junction.
    """
    covered = [c for c, (_, depth, _) in enumerate(columns) if depth]
    if not covered:
        return ""
    lo, hi = covered[0], covered[-1] + 1
    width = max(len("consensus"), *(len(n) for n in aligned)) + 2
    consensus = "".join(symbol or " " for symbol, _, _ in columns[lo:hi])
    lines = []
    if junction_col is not None and lo <= junction_col < hi:
        lines.append(" " * width + " " * (junction_col - lo) + "|")
    lines.append("consensus".ljust(width) + consensus)
    if EXPECTED in aligned:
        lines.append(EXPECTED.ljust(width) + aligned[EXPECTED][lo:hi])
    reads = sorted((n for n in aligned if n != EXPECTED), key=lambda n: _span(aligned[n]))
    for name in reads:
        row = aligned[name]
        a, b = _span(row)
        chars = []
        for c in range(lo, hi):
            if c < a or c > b:
                chars.append(" ")
            elif row[c] == columns[c][0]:
                chars.append(".")
            else:
                chars.append(row[c])
        lines.append(name.ljust(width) + "".join(chars))
    return "\n".join(lines) + "\n"


def align_junction(reads: dict[str, str], split_offsets: dict[str, int],
                   expected: tuple[str, int] | None) -> dict:
    """Align *reads* with MAFFT and summarise the junction.

    *split_offsets* gives, for split reads, the index of their first base right
    of the junction; *expected* is ``(sequence, left_length)`` or ``None``.
    """
    expected_seq, expected_offset = expected if expected else (None, None)
    aligned = run_mafft(reads, expected_seq)
    columns = consensus_columns(aligned)
    junction_col = junction_column(aligned, expected_offset, split_offsets)
    return {
        "aligned": aligned,
        "columns": columns,
        "junction_col": junction_col,
        "consensus": consensus_sequence(columns),
        "stats": junction_stats(aligned, columns, junction_col, expected_offset),
        "text": alignment_text(aligned, columns, junction_col),
    }
