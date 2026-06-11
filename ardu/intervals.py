"""
Genomic interval / region utilities.

Covers:
    - parsing a 4-column BED-like region file
    - computing the total span of a multi-region locus
    - per-locus custom plot-interval parsing
    - in-region coverage median helper
"""

from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# BED-like region file parsing
# ---------------------------------------------------------------------------

def get_regions(region_file: str) -> dict[str, list[str]]:
    """Parse a 4-column BED-like file into ``{locus_name: [region_str, ...]}``.

    File format (tab-delimited, no header)::

        chr1    10000    20000    LOCUS_A
        chr1    20500    21000    LOCUS_A
        chr2    50000    60000    LOCUS_B

    Lines starting with ``#`` or empty lines are silently ignored.
    Lines that do not contain exactly 4 tab-delimited fields are printed
    as warnings and skipped.

    Parameters
    ----------
    region_file:
        Path to the BED-like file.

    Returns
    -------
    dict[str, list[str]]
        Keys are locus names; values are lists of SAM-style region strings
        (``"chrom:start-end"``).
    """
    regions_dict: dict[str, list[str]] = {}
    with open(region_file) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.strip().split("\t")
            if len(fields) != 4:
                print(f"Ignoring invalid line: {fields}")
                continue
            region_str = f"{fields[0]}:{fields[1]}-{fields[2]}"
            regions_dict.setdefault(fields[3], []).append(region_str)
    return regions_dict


# ---------------------------------------------------------------------------
# Locus span
# ---------------------------------------------------------------------------

def get_totalspan(
    regions_dict: dict[str, list[str]],
    locus_name: str,
) -> tuple[str, int, int] | None:
    """Return ``(chrom, min_start, max_end)`` for *locus_name*.

    Uses the pre-parsed *regions_dict* rather than re-reading the file.

    Returns ``None`` if *locus_name* is not present in the dict.
    """
    region_strings = regions_dict.get(locus_name)
    if not region_strings:
        return None

    starts, ends, chroms = [], [], []
    for r in region_strings:
        chrom, coords = r.split(":")
        s, e = coords.split("-")
        chroms.append(chrom)
        starts.append(int(s))
        ends.append(int(e))

    return chroms[0], min(starts), max(ends)


# ---------------------------------------------------------------------------
# Custom plot interval file
# ---------------------------------------------------------------------------

def parse_plot_intervals(file_path: str) -> dict[str, str]:
    """Parse a tab-delimited file of per-locus custom plot intervals.

    File format::

        LOCUS_A    chr1:9000-22000
        LOCUS_B    chr2:48000-62000

    Returns
    -------
    dict[str, str]
        Mapping locus name → SAM-style region string.
    """
    plot_intervals: dict[str, str] = {}
    with open(file_path) as fh:
        for line in fh:
            locus, interval = line.strip().split("\t")
            plot_intervals[locus] = interval
    return plot_intervals


# ---------------------------------------------------------------------------
# In-region median
# ---------------------------------------------------------------------------

def region_median(
    coverage: dict[int, float],
    start: int,
    end: int,
) -> float:
    """Compute median coverage in ``[start, end)`` from a coverage dict.

    Returns ``0`` if no positions fall within the range.
    """
    vals = [coverage[pos] for pos in range(start, end) if pos in coverage]
    return float(np.median(vals)) if vals else 0.0
