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
            end=pos,
            quality_threshold=13,
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
                coverage_values: list[int] = []
                for region in region_list:
                    if not region:
                        continue
                    for pileupcolumn in samfile.pileup(
                        region=region,
                        min_base_quality=13,
                        min_mapping_quality=20,
                    ):
                        coverage_values.append(pileupcolumn.n)

                if coverage_values:
                    mean_cov = np.mean(coverage_values)
                    median_cov = np.median(coverage_values)
                    sd_cov = (
                        np.std(coverage_values, ddof=1)
                        if len(coverage_values) > 1
                        else 0.0
                    )
                    covbases = len(coverage_values)
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
