"""
Breakpoint detection methods.

Covers:
    - univariate ruptures (per-sample)
    - multivariate ruptures (pooled across samples)
    - rolling-average shift detection (wraps ardu.signal.detect_shifts)
"""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd
import ruptures as rpt

from ardu.signal import detect_shifts


# ---------------------------------------------------------------------------
# Univariate ruptures (per-sample)
# ---------------------------------------------------------------------------

def run_ruptures_univariate(
    d: pd.DataFrame,
    algo_name: str = "BottomUp",
    model: str = "l2",
    signal_col: str = "moving_average",
    n_bkps: int = 2,
    pen: int | None = None,
) -> list[dict]:
    """Run ruptures on a single sample's signal.

    Parameters
    ----------
    d:
        Per-base depth DataFrame with at least ``"pos"`` and *signal_col*.
    algo_name:
        ruptures algorithm class name (e.g. ``"BottomUp"``, ``"Pelt"``).
    model:
        ruptures cost model (e.g. ``"l2"``, ``"rbf"``).
    signal_col:
        Column of *d* to use as the signal.
    n_bkps:
        Expected number of breakpoints (used when *pen* is ``None``).
    pen:
        Penalty value passed to ``predict(pen=pen)``; overrides *n_bkps*
        when provided.

    Returns
    -------
    list[dict]
        Each dict has keys ``"row_index"`` (0-based position in *d*) and
        ``"breakpoint_number"`` (1-based counter).  The calling code is
        responsible for resolving the genomic position from ``d.loc[row_index]``.
    """
    algo_class = getattr(rpt, algo_name)
    signal = d[signal_col].to_numpy().reshape(-1, 1)
    algo = algo_class(model=model).fit(signal)
    result = algo.predict(pen=pen) if pen is not None else algo.predict(n_bkps=n_bkps)

    breakpoints = []
    for line_num, b in enumerate(result[:-1], start=1):
        breakpoints.append({"row_index": b, "breakpoint_number": line_num})
    return breakpoints


# ---------------------------------------------------------------------------
# Multivariate ruptures (pooled across samples)
# ---------------------------------------------------------------------------

def run_ruptures_multivariate(
    sample_list: list[tuple[str, pd.DataFrame, str]],
    signal_col: str = "moving_average",
    algo_name: str = "BottomUp",
    model: str = "l2",
    n_bkps: int = 2,
    pen: int | None = None,
) -> tuple[np.ndarray, list[dict]]:
    """Run multivariate ruptures on a pooled set of samples.

    Samples are aligned to a shared position array (union of all positions)
    and interpolated to fill gaps before stacking into a signal matrix.

    Parameters
    ----------
    sample_list:
        List of ``(bam_name, depth_df, span_str)`` tuples, as stored in
        ``pooled_depth``.
    signal_col:
        Column name used as the signal from each sample's DataFrame.
    algo_name, model, n_bkps, pen:
        Passed through to ruptures.

    Returns
    -------
    ref_pos : np.ndarray
        The shared genomic position array (union of all sample positions).
    breakpoints : list[dict]
        Each dict has keys ``"position"`` (genomic bp), ``"breakpoint_number"``,
        and ``"method": "ruptures_pooled_multivariate"``.
    """
    if not sample_list:
        return np.array([]), []

    # Build union position array
    all_pos: set[int] = set()
    for _bam_name, d, _span in sample_list:
        all_pos.update(d["pos"].values)
    ref_pos = np.array(sorted(all_pos))

    def _dedup_and_align(d: pd.DataFrame, col: str) -> np.ndarray:
        s = d.groupby("pos", sort=True)[col].mean()
        s = s.reindex(ref_pos).interpolate(limit_direction="both")
        return s.values

    aligned = [_dedup_and_align(d, signal_col) for _bam_name, d, _span in sample_list]
    signal_matrix = np.column_stack(aligned)

    # Replace any remaining NaNs with column means
    col_means = np.nanmean(signal_matrix, axis=0)
    nan_idx = np.where(np.isnan(signal_matrix))
    signal_matrix[nan_idx] = np.take(col_means, nan_idx[1])

    algo_class = getattr(rpt, algo_name)
    algo = algo_class(model=model).fit(signal_matrix)
    result = algo.predict(pen=pen) if pen is not None else algo.predict(n_bkps=n_bkps)

    breakpoints = []
    for line_num, b in enumerate(result[:-1], start=1):
        try:
            bp_pos = int(ref_pos[b])
        except (IndexError, ValueError, TypeError):
            print(f"Warning: invalid position at row {b} in pooled ruptures, skipping.")
            continue
        breakpoints.append(
            {
                "position": bp_pos,
                "breakpoint_number": line_num,
                "method": "ruptures_pooled_multivariate",
            }
        )
    return ref_pos, breakpoints


# ---------------------------------------------------------------------------
# Rolling-average shift detection
# ---------------------------------------------------------------------------

def run_rollingaverage(
    d: pd.DataFrame,
    threshold: float = 0.5,
    window_size: int = 1000,
    passes: int = 1,
) -> list[tuple[int, int]]:
    """Detect breakpoint regions by rolling-average shifts.

    Wraps :func:`ardu.signal.detect_shifts` and then merges consecutive
    flagged positions into contiguous regions.

    Parameters
    ----------
    d:
        Per-base depth DataFrame; must contain ``"pos"``, ``"norm"``, and
        ``"moving_average"``.
    threshold:
        Absolute signal–MA difference that triggers a flag.
    window_size:
        Window width for extra smoothing inside :func:`detect_shifts`.
    passes:
        Number of smoothing passes.

    Returns
    -------
    list[tuple[int, int]]
        List of ``(region_start, region_end)`` genomic coordinate pairs.
    """
    significant_shifts = detect_shifts(d, threshold, window_size, passes)
    regions_bkp: list[tuple[int, int]] = []
    consecutive_positions: list[int] = []

    for _bidx, brow in significant_shifts.iterrows():
        if not consecutive_positions or brow["pos"] - consecutive_positions[-1] <= 1:
            consecutive_positions.append(brow["pos"])
        else:
            regions_bkp.append((consecutive_positions[0], consecutive_positions[-1]))
            consecutive_positions = [brow["pos"]]

    if consecutive_positions:
        regions_bkp.append((consecutive_positions[0], consecutive_positions[-1]))

    return regions_bkp
