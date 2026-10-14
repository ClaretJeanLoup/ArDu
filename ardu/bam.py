"""
BAM / pysam I/O utilities.

Covers:
    - per-position soft-clip counting
    - per-nucleotide read counts at a single position
    - depth-of-coverage statistics across a set of named loci
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm


# ---------------------------------------------------------------------------
# Shared depth computation
# ---------------------------------------------------------------------------

MIN_BASE_QUALITY = 13
MIN_MAPPING_QUALITY = 20


def passes_read_filters(read: pysam.AlignedSegment) -> bool:
    """Read-level filter shared by every depth computation in ArDu."""
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_qcfail
        or read.is_duplicate
    ) and read.mapping_quality >= MIN_MAPPING_QUALITY


def parse_region(region: str) -> tuple[str, int, int]:
    """Parse a samtools-style ``chrom:start-end`` (1-based, inclusive) region.

    Returns ``(chrom, start, end)`` as 0-based half-open coordinates.
    """
    chrom, _, span = region.rpartition(":")
    if not chrom:
        raise ValueError(f"Region '{region}' has no start-end range")
    start_str, end_str = span.replace(",", "").split("-")
    return chrom, int(start_str) - 1, int(end_str)


def depth_array(
    samfile: pysam.AlignmentFile,
    chrom: str,
    start: int,
    end: int,
) -> np.ndarray:
    """Per-base depth over ``[start, end)``, 0-based, with zeros where uncovered.

    Bases below ``MIN_BASE_QUALITY`` and reads failing
    :func:`passes_read_filters` are not counted.
    """
    start = max(0, start)
    end = min(end, samfile.get_reference_length(chrom))
    if end <= start:
        return np.zeros(0, dtype=np.int64)
    acgt = samfile.count_coverage(
        chrom,
        start=start,
        stop=end,
        quality_threshold=MIN_BASE_QUALITY,
        read_callback=passes_read_filters,
    )
    return np.sum(np.array(acgt), axis=0)


def depth_over_region(bam_file: str, region: str) -> list[tuple[int, int]]:
    """``(pos, depth)`` for every base of a ``chrom:start-end`` region.

    Positions are 0-based; uncovered bases are reported with depth 0.
    """
    chrom, start, end = parse_region(region)
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        depth = depth_array(samfile, chrom, start, end)
    start = max(0, start)
    return list(zip(range(start, start + len(depth)), depth.tolist()))


# ---------------------------------------------------------------------------
# Soft-clip detection
# ---------------------------------------------------------------------------

def detect_softclip_counts(
    bam_file: str,
    region: str,
    min_clip: int = 5,
) -> pd.DataFrame:
    """Count clipped reads per position across *region*.

    Parameters
    ----------
    bam_file:
        Path to an indexed BAM file.
    region:
        SAM-style region string, e.g. ``"chr1:100000-200000"``.
    min_clip:
        Minimum number of soft-clipped bases to count a read.

    Returns
    -------
    pd.DataFrame
        Columns: ``chromosome``, ``position``, ``softclip_count``.
        Sorted by position.
    """
    softclip_counts: dict[int, int] = {}
    chrom = region.split(":")[0]

    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        for read in samfile.fetch(region=region):
            if read.cigartuples is None:
                continue
            # Left-side clip
            if (
                read.cigartuples[0][0] == 4
                and read.cigartuples[0][1] >= min_clip
            ):
                pos = read.reference_start
                softclip_counts[pos] = softclip_counts.get(pos, 0) + 1
            # Right-side clip
            if (
                read.cigartuples[-1][0] == 4
                and read.cigartuples[-1][1] >= min_clip
            ):
                pos = read.reference_end
                softclip_counts[pos] = softclip_counts.get(pos, 0) + 1

    return pd.DataFrame(
        [
            {"chromosome": chrom, "position": pos, "softclip_count": count}
            for pos, count in softclip_counts.items()
        ]
    ).sort_values("position").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Nucleotide counts at a single position
# ---------------------------------------------------------------------------

def calculate_nucleotide_counts(
    bam_file: str,
    chromosome: str,
    position: int | str,
) -> tuple[dict[str, int], int]:
    """Return per-nucleotide read counts and total depth at *position*.

    Parameters
    ----------
    bam_file:
        Path to an indexed BAM file.
    chromosome:
        Chromosome / contig name.
    position:
        1-based genomic position (int or string).

    Returns
    -------
    counts : dict[str, int]
        Keys ``"A"``, ``"C"``, ``"G"``, ``"T"`` — base counts at *position*.
    total_depth : int
        Sum of all four nucleotide counts.
    """
    pos = int(position)
    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        a, c, g, t = samfile.count_coverage(
            chromosome,
            start=pos - 1,
            stop=pos,
            quality_threshold=MIN_BASE_QUALITY,
            read_callback=passes_read_filters,
        )
    counts = {
        "A": int(a[0]) if len(a) > 0 else 0,
        "C": int(c[0]) if len(c) > 0 else 0,
        "G": int(g[0]) if len(g) > 0 else 0,
        "T": int(t[0]) if len(t) > 0 else 0,
    }
    return counts, sum(counts.values())


# ---------------------------------------------------------------------------
# Coverage statistics across loci
# ---------------------------------------------------------------------------

def calculate_coverage_stats(
    regions: dict[str, list[str]],
    bam_file: str,
) -> pd.DataFrame:
    """Compute mean, median, and sd of depth of coverage for each named locus.

    Parameters
    ----------
    regions:
        Dict mapping locus name → list of SAM region strings, as returned by
        :func:`ardu.intervals.get_regions`.
    bam_file:
        Path to an indexed BAM file.

    Returns
    -------
    pd.DataFrame
        Columns: ``locus``, ``mean``, ``median``, ``sd``, ``coveredbases``.
        Numeric values where coverage data exists; the string ``"NA"``
        otherwise.
    """
    data_list: list[dict] = []
    try:
        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            for name, region_list in tqdm(
                regions.items(),
                total=len(regions),
                desc="Processing targets",
            ):
                per_region = []
                for region in region_list:
                    if not region:
                        continue
                    chrom, start, end = parse_region(region)
                    per_region.append(depth_array(samfile, chrom, start, end))
                coverage_values = (
                    np.concatenate(per_region) if per_region else np.zeros(0)
                )

                if coverage_values.size:
                    mean_cov = np.mean(coverage_values)
                    median_cov = np.median(coverage_values)
                    sd_cov = (
                        np.std(coverage_values, ddof=1)
                        if coverage_values.size > 1
                        else 0.0
                    )
                    covbases = int(np.count_nonzero(coverage_values))
                else:
                    mean_cov = median_cov = sd_cov = covbases = "NA"

                data_list.append(
                    {
                        "locus": name,
                        "mean": mean_cov,
                        "median": median_cov,
                        "sd": sd_cov,
                        "coveredbases": covbases,
                    }
                )
    except Exception as exc:
        print(f"Error processing BAM file {bam_file}: {exc}")
        return pd.DataFrame()

    return pd.DataFrame(data_list)
