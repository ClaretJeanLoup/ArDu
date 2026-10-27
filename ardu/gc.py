"""
GC content of loci, and random normalisation windows matched on GC.

Covers:
    - GC of a set of regions, padded on each side to the scale of a fragment
    - the GC report comparing each locus with the normalisation locus
    - drawing random windows away from the loci, with a GC close to theirs
"""

from __future__ import annotations

import bisect

import numpy as np
import pysam

from ardu.bam import parse_region


# fragments covering a base mostly overlap its close neighbours: across 46 Anopheles bams, 100 bp of padding
# predicted the depth of a reference locus from the GC curve best among 0, 50, 100, 150 and 250 bp
GC_PAD = 100
MATCH_TOLERANCE = 0.02
MIN_CONTIG = 1_000_000
MAX_N_FRACTION = 0.1


# ---------------------------------------------------------------------------
# GC of regions
# ---------------------------------------------------------------------------

def _merge(intervals: list[tuple[int, int]]) -> list[list[int]]:
    merged: list[list[int]] = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def region_gc(
    fasta: pysam.FastaFile,
    regions: list[tuple[str, int, int]],
    pad: int = GC_PAD,
) -> tuple[float, int]:
    """GC fraction of *regions* (0-based half-open), each padded by *pad* bp.

    Overlapping padded regions are counted once. Returns ``(gc, n)`` where *n*
    is the number of A, C, G or T bases counted; *gc* is NaN when *n* is 0.
    """
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in regions:
        length = fasta.get_reference_length(chrom)
        by_chrom.setdefault(chrom, []).append((max(0, start - pad), min(length, end + pad)))
    gc = acgt = 0
    for chrom, intervals in by_chrom.items():
        for s, e in _merge(intervals):
            seq = fasta.fetch(chrom, s, e).upper()
            g = seq.count("G") + seq.count("C")
            gc += g
            acgt += g + seq.count("A") + seq.count("T")
    return (gc / acgt if acgt else float("nan")), acgt


def loci_gc(
    fasta_file: str,
    regions_dict: dict[str, list[str]],
    pad: int = GC_PAD,
) -> dict[str, tuple[float, int]]:
    """``{locus: (gc, n_bases)}`` for every locus of *regions_dict*.

    Loci on contigs missing from the fasta get ``(nan, 0)``.
    """
    out = {}
    with pysam.FastaFile(fasta_file) as fasta:
        for name, region_list in regions_dict.items():
            regions = [parse_region(r) for r in region_list if r]
            if any(chrom not in fasta.references for chrom, _, _ in regions):
                out[name] = (float("nan"), 0)
            else:
                out[name] = region_gc(fasta, regions, pad)
    return out


def write_gc_report(
    path: str,
    gc_by_locus: dict[str, tuple[float, int]],
    norm: str,
    warn_points: float,
    pad: int = GC_PAD,
) -> list[str]:
    """Write the GC of each locus and its difference with *norm*, in GC points.

    Returns the loci whose difference exceeds *warn_points*.
    """
    norm_gc = gc_by_locus[norm][0]
    flagged = []
    with open(path, "w") as fh:
        fh.write("# ArDu GC report\n")
        fh.write(f"# gc: GC of the locus regions padded by {pad} bp on each side (A, C, G, T bases only).\n")
        fh.write(f"# diff: locus GC minus normalisation locus ({norm}) GC, in percentage points; "
                 f"warn when above {warn_points:g}.\n")
        fh.write("locus\tbases\tgc\tdiff\twarning\n")
        for name, (gc, n) in gc_by_locus.items():
            diff = (gc - norm_gc) * 100
            warn = name != norm and not np.isnan(diff) and abs(diff) > warn_points
            if warn:
                flagged.append(name)
            gc_str = "NA" if np.isnan(gc) else f"{gc:.4f}"
            diff_str = "NA" if np.isnan(diff) else f"{diff:.1f}"
            fh.write(f"{name}\t{n}\t{gc_str}\t{diff_str}\t{'yes' if warn else ''}\n")
    return flagged


# ---------------------------------------------------------------------------
# Random normalisation windows
# ---------------------------------------------------------------------------

def read_bed_intervals(path: str) -> list[tuple[str, int, int]]:
    """First three columns of a BED file as 0-based half-open intervals."""
    out = []
    with open(path) as fh:
        for raw in fh:
            if not raw.strip() or raw.startswith(("#", "track", "browser")):
                continue
            f = raw.split("\t")
            out.append((f[0], int(f[1]), int(f[2])))
    return out


def draw_norm_windows(
    fasta_file: str,
    n: int,
    size: int,
    exclude: list[tuple[str, int, int]],
    target_gc: float | None = None,
    chroms: list[str] | None = None,
    seed: int = 1,
    pad: int = GC_PAD,
    tolerance: float = MATCH_TOLERANCE,
) -> list[tuple[str, int, int]]:
    """Draw up to *n* non-overlapping windows of *size* bp, sorted by position.

    Windows are drawn uniformly over *chroms* (default: contigs of at least
    ``MIN_CONTIG`` bp), outside *exclude*, with at most ``MAX_N_FRACTION`` of
    N. With *target_gc*, only windows whose padded GC is within *tolerance*
    of it are kept. Stops after ``200 * n`` candidates, so fewer than *n*
    windows come back when few parts of the genome match.
    """
    rng = np.random.default_rng(seed)
    with pysam.FastaFile(fasta_file) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        if chroms:
            missing = [c for c in chroms if c not in lengths]
            if missing:
                raise ValueError(f"Contigs not in the fasta: {', '.join(missing)}")
            contigs = [c for c in chroms if lengths[c] > size]
        else:
            contigs = [c for c, l in lengths.items() if l >= MIN_CONTIG]
        if not contigs:
            raise ValueError("No contig to draw normalisation windows from.")
        # windows sit on a grid of their own size, so two windows never overlap
        n_slots = np.array([lengths[c] // size for c in contigs])
        weights = n_slots / n_slots.sum()

        blocked: dict[str, tuple[list[int], list[int]]] = {}
        for chrom, s, e in sorted(exclude):
            starts, ends = blocked.setdefault(chrom, ([], []))
            if starts and s <= ends[-1]:
                ends[-1] = max(ends[-1], e)
            else:
                starts.append(s)
                ends.append(e)

        def overlaps_excluded(chrom: str, s: int, e: int) -> bool:
            if chrom not in blocked:
                return False
            starts, ends = blocked[chrom]
            i = bisect.bisect_right(starts, e - 1) - 1
            return i >= 0 and ends[i] > s

        chosen: set[tuple[int, int]] = set()
        tried: set[tuple[int, int]] = set()
        for _ in range(200):
            ci = rng.choice(len(contigs), size=n, p=weights)
            slots = (rng.random(n) * n_slots[ci]).astype(int)
            for c, k in zip(ci.tolist(), slots.tolist()):
                if (c, k) in tried:
                    continue
                tried.add((c, k))
                chrom, s = contigs[c], k * size
                e = s + size
                if overlaps_excluded(chrom, s, e):
                    continue
                seq = fasta.fetch(chrom, s, e).upper()
                if seq.count("N") > MAX_N_FRACTION * size:
                    continue
                if target_gc is not None:
                    gc, _ = region_gc(fasta, [(chrom, s, e)], pad)
                    if np.isnan(gc) or abs(gc - target_gc) > tolerance:
                        continue
                chosen.add((c, k))
                if len(chosen) == n:
                    break
            if len(chosen) == n:
                break
    return sorted((contigs[c], k * size, k * size + size) for c, k in chosen)


def write_windows_bed(path: str, windows: list[tuple[str, int, int]], name: str) -> None:
    """Write *windows* in the 4-column region format (1-based starts), all named *name*."""
    with open(path, "w") as fh:
        for chrom, s, e in windows:
            fh.write(f"{chrom}\t{s + 1}\t{e}\t{name}\n")
