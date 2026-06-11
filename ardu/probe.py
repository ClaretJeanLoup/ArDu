"""
Probe-based interval utilities.

Covers:
    - adaptive interval expansion (expand_interval)
    - per-sample probe-profile building (build_probe_profile)
    - probe-profile clustering to group structurally similar samples
      (cluster_probe_profiles)
    - profile serialisation helpers
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
import pysam


# ---------------------------------------------------------------------------
# Adaptive interval expansion
# ---------------------------------------------------------------------------

def expand_interval(
    bam_file: str,
    chrom: str,
    start: int,
    end: int,
    target_cov: float,
    reference_cov: float,
    probe_size: int = 100,
    probe_spacing: int | None = None,
    probe_number: int = 20,
    stop_drops: int = 7,
    max_extension: int | None = None,
    verbose: bool = True,
    drop_threshold: float = 5,
    null_threshold: int = 3,
    min_extension_factor: float = 2,
    symmetric: bool = False,
    max_spacing: int | None = None,
) -> str:
    """Adaptively expand a genomic interval outward until coverage drops.

    Probes are placed at regularly-spaced positions to the left and right of
    ``[start, end]``.  Expansion stops when *stop_drops* consecutive probes
    show coverage below *reference_cov*.

    Parameters
    ----------
    bam_file:
        Path to an indexed BAM.
    chrom:
        Chromosome name.
    start, end:
        0-based half-open coordinates of the target interval.
    target_cov:
        Mean/median coverage inside the target (used for labelling / debug).
    reference_cov:
        Background coverage; a probe below this value counts as a "drop".
    probe_size:
        Length of each probe window in bp.
    probe_spacing:
        Distance between probe starts; defaults to *probe_size* if ``None``.
    probe_number:
        Number of probes evaluated per round.
    stop_drops:
        Minimum run of consecutive below-threshold probes that signals
        "coverage has returned to baseline — stop expanding".
    max_extension:
        Hard limit on how far each side can extend from the original
        boundary (bp).  ``None`` means 10× the target length.
    verbose:
        Print debug lines.
    drop_threshold:
        Not used directly; kept for API compatibility.
    null_threshold:
        Number of zero-coverage probes that triggers a coarser jump.
    min_extension_factor:
        The final interval must be at least this multiple of the original
        target length; if not, it is expanded symmetrically to meet the
        requirement.
    symmetric:
        Force equal extension on both sides.
    max_spacing:
        Cap on the adaptive probe spacing growth.

    Returns
    -------
    str
        SAM-style region string ``"chrom:left-right"`` for the expanded
        interval.
    """
    if probe_spacing is None:
        probe_spacing = probe_size
    if reference_cov is None:
        raise ValueError("reference_cov must be provided.")

    with pysam.AlignmentFile(bam_file, "rb") as samfile:
        original_left, original_right = start, end
        target_len = original_right - original_left
        left_max_ext = (
            original_left - max_extension
            if max_extension
            else max(0, original_left - target_len * 10)
        )
        right_max_ext = (
            original_right + max_extension
            if max_extension
            else original_right + target_len * 10
        )
        fetch_start = max(0, left_max_ext)
        fetch_end = right_max_ext
        if verbose:
            print(f"[DEBUG] Fetching coverage for {chrom}:{fetch_start}-{fetch_end} ...")
        cov_tuple = samfile.count_coverage(
            chrom, start=fetch_start, end=fetch_end, quality_threshold=13
        )
        coverage_array = np.sum(np.array(cov_tuple), axis=0)

    if verbose:
        print(f"[DEBUG] Coverage fetched, length={len(coverage_array)} bases")

    offset = fetch_start
    left_boundary = original_left
    right_boundary = original_right
    left_done = right_done = False
    rounds = rounds_no_drop_left = rounds_no_drop_right = 0

    def _has_consecutive(probe_covs: np.ndarray, threshold: float, min_run: int) -> bool:
        run = 0
        for c in probe_covs:
            if c < threshold:
                run += 1
                if run >= min_run:
                    return True
            else:
                run = 0
        return False

    while not (left_done and right_done):
        rounds += 1
        left_spacing = probe_spacing * (1 + rounds_no_drop_left)
        right_spacing = probe_spacing * (1 + rounds_no_drop_right)
        if max_spacing:
            left_spacing = min(left_spacing, max_spacing)
            right_spacing = min(right_spacing, max_spacing)

        left_starts = [
            max(0, int(left_boundary - (i + 1) * left_spacing - offset))
            for i in range(probe_number)
        ]
        left_ends = [s + probe_size for s in left_starts]
        right_starts = [
            max(0, int(right_boundary + i * right_spacing - offset))
            for i in range(probe_number)
        ]
        right_ends = [s + probe_size for s in right_starts]

        left_covs = np.array(
            [coverage_array[s:e].mean() if e > s else 0 for s, e in zip(left_starts, left_ends)]
        )
        right_covs = np.array(
            [coverage_array[s:e].mean() if e > s else 0 for s, e in zip(right_starts, right_ends)]
        )

        left_drops = np.sum((left_covs < reference_cov) & (left_covs > 0))
        right_drops = np.sum((right_covs < reference_cov) & (right_covs > 0))
        left_nulls = np.sum(left_covs == 0)
        right_nulls = np.sum(right_covs == 0)

        stop_left = _has_consecutive(left_covs, target_cov, stop_drops)
        stop_right = _has_consecutive(right_covs, target_cov, stop_drops)

        if verbose:
            print(
                f"[DEBUG] Round {rounds}: L={left_boundary}, R={right_boundary}\n"
                f"        left_drops={left_drops}, right_drops={right_drops}, "
                f"left_nulls={left_nulls}, right_nulls={right_nulls}, "
                f"expand_left={not stop_left}, expand_right={not stop_right}"
            )

        if not stop_left:
            if left_nulls >= null_threshold:
                left_boundary = max(0, left_boundary - 2 * target_len)
            else:
                expand_amount = max(1, int(probe_number * left_spacing))
                left_boundary = max(0, left_boundary - expand_amount)
                rounds_no_drop_left = (
                    rounds_no_drop_left + 1 if left_drops == 0 else 0
                )
        else:
            rounds_no_drop_left = 0

        if not stop_right:
            if right_nulls >= null_threshold:
                right_boundary += 2 * target_len
            else:
                expand_amount = max(1, int(probe_number * right_spacing))
                right_boundary += expand_amount
                rounds_no_drop_right = (
                    rounds_no_drop_right + 1 if right_drops == 0 else 0
                )
        else:
            rounds_no_drop_right = 0

        if max_extension:
            left_boundary = max(original_left - max_extension, left_boundary)
            right_boundary = min(original_right + max_extension, right_boundary)

        left_done = stop_left
        right_done = stop_right

        if rounds > 1_000_000:
            if verbose:
                print("[WARNING] Maximum rounds exceeded, stopping.")
            break

    if symmetric:
        left_ext = original_left - left_boundary
        right_ext = right_boundary - original_right
        max_ext = max(left_ext, right_ext)
        left_boundary = original_left - max_ext
        right_boundary = original_right + max_ext

    min_extension = (original_right - original_left) * min_extension_factor
    if (right_boundary - left_boundary) < min_extension:
        mid = (original_left + original_right) // 2
        half = int(min_extension // 2)
        left_boundary = max(0, mid - half)
        right_boundary = mid + half
        if max_extension:
            left_boundary = max(original_left - max_extension, left_boundary)
            right_boundary = min(original_right + max_extension, right_boundary)

    if verbose:
        print(
            f"[RESULT] Final interval: {chrom}:{left_boundary}-{right_boundary} | Rounds={rounds}"
        )
    return f"{chrom}:{left_boundary}-{right_boundary}"


# ---------------------------------------------------------------------------
# Probe-profile building
# ---------------------------------------------------------------------------

def build_probe_profile(
    bam_file: str,
    chrom: str,
    start: int,
    end: int,
    ref_median: float,
    probe_size: int = 500,
    probe_spacing: int = 1000,
    n_probes_each_side: int = 20,
) -> tuple[np.ndarray, np.ndarray]:
    """Launch evenly-spaced probes symmetrically around ``[start, end]``.

    Parameters
    ----------
    bam_file:
        Path to an indexed BAM.
    chrom, start, end:
        Target locus coordinates (0-based half-open).
    ref_median:
        Normalisation coverage; probe ratios are expressed as
        ``probe_cov / ref_median``.
    probe_size:
        Width of each probe window in bp.
    probe_spacing:
        Centre-to-centre distance between consecutive probes in bp.
    n_probes_each_side:
        Number of probes placed on each side of the target (and across
        the target centre), giving ``2 * n_probes_each_side + 1`` total.

    Returns
    -------
    positions : np.ndarray, shape ``(n_probes,)``
        Genomic centre of each probe.
    ratios : np.ndarray, shape ``(n_probes,)``
        ``probe_cov / ref_median`` for each probe.
    """
    fetch_start = max(
        0, start - n_probes_each_side * probe_spacing - probe_size
    )
    fetch_end = end + n_probes_each_side * probe_spacing + probe_size

    with pysam.AlignmentFile(bam_file, "rb") as sam:
        chrom_len = dict(zip(sam.references, sam.lengths)).get(chrom, fetch_end)
        fetch_end = min(fetch_end, chrom_len)
        raw = sam.count_coverage(
            chrom, start=fetch_start, end=fetch_end, quality_threshold=13
        )
    coverage_array = np.sum(np.array(raw), axis=0)
    offset = fetch_start
    target_centre = (start + end) // 2

    positions: list[int] = []
    ratios: list[float] = []

    for i in range(-n_probes_each_side, n_probes_each_side + 1):
        probe_centre = target_centre + i * probe_spacing
        probe_s = max(0, probe_centre - probe_size // 2 - offset)
        probe_e = min(probe_s + probe_size, len(coverage_array))
        probe_cov = (
            coverage_array[probe_s:probe_e].mean() if probe_e > probe_s else 0.0
        )
        positions.append(probe_centre)
        ratios.append(probe_cov / ref_median if ref_median > 0 else 0.0)

    return np.array(positions), np.array(ratios)


# ---------------------------------------------------------------------------
# Probe-profile clustering
# ---------------------------------------------------------------------------

def cluster_probe_profiles(
    profiles: dict[str, np.ndarray],
    threshold: float = 1.4,
    distance_threshold: float = 0.15,
    method: str = "hamming",
) -> dict[str, int]:
    """Group samples by structural similarity of their probe profiles.

    Each profile is first binarised (probe > *threshold* → 1, else 0),
    then samples are clustered by Hamming distance or exact string match.

    Parameters
    ----------
    profiles:
        ``{bam_name: ratio_array}`` as returned by
        :func:`build_probe_profile`.
    threshold:
        Coverage ratio above which a probe position is considered
        "elevated" (copy-number gain).
    distance_threshold:
        Maximum Hamming distance for two samples to be placed in the
        same cluster.  Only used when ``method="hamming"``.
    method:
        ``"hamming"`` (default, uses scikit-learn) or ``"string"``
        (exact binary string match, no external dependencies).

    Returns
    -------
    dict[str, int]
        ``{bam_name: group_id}``
    """
    names = list(profiles.keys())
    if not names:
        return {}

    binary = np.vstack([(profiles[n] > threshold).astype(int) for n in names])

    if method == "string":
        profile_strings = ["".join(map(str, row)) for row in binary]
        unique = {s: i for i, s in enumerate(dict.fromkeys(profile_strings))}
        return {name: unique[ps] for name, ps in zip(names, profile_strings)}

    # Hamming-distance agglomerative clustering
    try:
        from sklearn.cluster import AgglomerativeClustering
    except ImportError as exc:
        raise ImportError(
            "scikit-learn is required for probe clustering. "
            "Install with: pip install scikit-learn"
        ) from exc

    if len(names) == 1:
        return {names[0]: 0}

    clustering = AgglomerativeClustering(
        n_clusters=None,
        distance_threshold=distance_threshold,
        metric="hamming",
        linkage="average",
        compute_full_tree=True,
    )
    labels_array = clustering.fit_predict(binary)
    return {name: int(lbl) for name, lbl in zip(names, labels_array)}


# ---------------------------------------------------------------------------
# Serialisation helpers
# ---------------------------------------------------------------------------

def profile_to_string(ratio_array: np.ndarray, threshold: float = 1.4) -> str:
    """Convert a ratio array to a human-readable binary string."""
    return "".join("1" if v > threshold else "0" for v in ratio_array)


def write_probe_profiles_tsv(
    profiles: dict[str, np.ndarray],
    positions: np.ndarray,
    labels: dict[str, int],
    threshold: float,
    output_path: str,
) -> None:
    """Write probe profiles and cluster labels to a TSV for inspection.

    Parameters
    ----------
    profiles:
        ``{bam_name: ratio_array}``
    positions:
        Genomic centres of each probe (same length as each ratio array).
    labels:
        ``{bam_name: group_id}`` as returned by
        :func:`cluster_probe_profiles`.
    threshold:
        Binarisation threshold used for the profile string column.
    output_path:
        Destination TSV file path.
    """
    rows = []
    for name, ratios in profiles.items():
        row: dict = {
            "sample": name,
            "group": labels.get(name, -1),
            "profile_string": profile_to_string(ratios, threshold),
        }
        for pos, ratio in zip(positions, ratios):
            row[f"pos_{pos}"] = round(float(ratio), 3)
        rows.append(row)
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)
