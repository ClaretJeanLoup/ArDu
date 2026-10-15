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

PENALTY_BIN = 500
REFINE_BINS = 5


def binned_signal(d: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Mean normalised depth per bin and the row of *d* where each bin starts.

    Bins are those of ``--plot-bin`` when *d* was binned, else
    :data:`PENALTY_BIN` bp.
    """
    keys = d["genomic_bins"] if "genomic_bins" in d.columns else d["pos"] // PENALTY_BIN
    means = d["norm"].groupby(keys, observed=True, sort=True).mean().to_numpy()
    starts = d.index.to_series(index=d.index).groupby(keys, observed=True, sort=True).min()
    return means, starts.to_numpy()


def scaled_penalty(signal: np.ndarray, pen: float) -> float:
    """*pen* times sigma squared times log(n), sigma being the noise of
    *signal* estimated from the differences between consecutive values.
    """
    diffs = np.diff(signal)
    sigma = 1.4826 * np.median(np.abs(diffs - np.median(diffs))) / np.sqrt(2)
    return pen * sigma ** 2 * np.log(len(signal))


def refine_breakpoint(norm: np.ndarray, start: int, end: int, min_size: int = 50) -> int:
    """Row in ``[start, end)`` that best splits *norm* there into two segments
    of different mean (least squares), at least *min_size* rows each.
    """
    x = norm[start:end]
    n = len(x)
    if n < 2 * min_size:
        return start + n // 2
    csum, csum2 = np.cumsum(x), np.cumsum(x ** 2)
    k = np.arange(min_size, n - min_size + 1)
    left = csum2[k - 1] - csum[k - 1] ** 2 / k
    right = (csum2[-1] - csum2[k - 1]) - (csum[-1] - csum[k - 1]) ** 2 / (n - k)
    return start + int(k[np.argmin(left + right)])


def run_ruptures_univariate(
    d: pd.DataFrame,
    algo_name: str = "BottomUp",
    model: str = "l2",
    signal_col: str = "moving_average",
    n_bkps: int = 2,
    pen: float | None = None,
) -> list[dict]:
    """Run ruptures on a single sample's signal.

    With *n_bkps*, ruptures looks for that many breakpoints in the per-base
    *signal_col*. With *pen*, the number of breakpoints is free: ruptures runs
    on the mean normalised depth per bin (see :func:`binned_signal`), and a
    breakpoint is kept when it lowers the cost by more than
    ``pen * sigma^2 * log(n)`` (:func:`scaled_penalty`). The penalty is then
    independent of depth, noise and region length; about 10 suits depth data,
    lower values find smaller steps. Each breakpoint is then moved to the best
    split of the per-base normalised depth within :data:`REFINE_BINS` bins on
    either side (:func:`refine_breakpoint`).

    Parameters
    ----------
    d:
        Per-base depth DataFrame with at least ``"pos"``, ``"norm"`` and *signal_col*.
    algo_name:
        ruptures algorithm class name (e.g. ``"BottomUp"``, ``"Pelt"``).
    model:
        ruptures cost model (e.g. ``"l2"``, ``"rbf"``).
    signal_col:
        Column of *d* to use as the signal (without *pen*).
    n_bkps:
        Expected number of breakpoints (used when *pen* is ``None``).
    pen:
        Penalty multiplier; overrides *n_bkps* when provided.

    Returns
    -------
    list[dict]
        Each dict has keys ``"row_index"`` (0-based position in *d*) and
        ``"breakpoint_number"`` (1-based counter).  The calling code is
        responsible for resolving the genomic position from ``d.loc[row_index]``.
    """
    algo_class = getattr(rpt, algo_name)
    if pen is None:
        signal = d[signal_col].to_numpy().reshape(-1, 1)
        result = algo_class(model=model).fit(signal).predict(n_bkps=n_bkps)
        rows = result[:-1]
    else:
        means, starts = binned_signal(d)
        if len(means) < 4:
            return []
        result = algo_class(model=model, min_size=2).fit(means.reshape(-1, 1)).predict(
            pen=scaled_penalty(means, pen)
        )
        norm = d["norm"].to_numpy()
        bounds = list(starts) + [len(d)]
        rows = [refine_breakpoint(norm, int(bounds[max(0, i - REFINE_BINS)]),
                                  int(bounds[min(len(bounds) - 1, i + REFINE_BINS)]))
                for i in result[:-1]]

    return [{"row_index": b, "breakpoint_number": line_num}
            for line_num, b in enumerate(rows, start=1)]


# ---------------------------------------------------------------------------
# Multivariate ruptures (pooled across samples)
# ---------------------------------------------------------------------------

def _pooled_penalty_breakpoints(
    sample_list: list[tuple[str, pd.DataFrame, str]],
    algo_name: str, model: str, pen: float,
) -> tuple[np.ndarray, list[dict]]:
    """Pooled breakpoints with a penalty instead of a set number.

    Each sample's mean normalised depth per :data:`PENALTY_BIN` bp bin (bins
    shared by all samples) is divided by its noise (see :func:`scaled_penalty`),
    and a breakpoint is kept when it lowers the summed cost by more than
    ``pen * n_samples * log(n_bins)``, the single-sample penalty summed over
    samples. Each breakpoint is then refined in every sample as in
    :func:`run_ruptures_univariate` and placed at the median position.
    """
    columns = []
    for _bam_name, d, _span in sample_list:
        pos_col = "base_pos" if "base_pos" in d.columns else "pos"
        means = d["norm"].groupby(d[pos_col] // PENALTY_BIN).mean()
        diffs = np.diff(means.to_numpy())
        sigma = 1.4826 * np.median(np.abs(diffs - np.median(diffs))) / np.sqrt(2)
        columns.append(means / (sigma or 1.0))
    matrix = pd.concat(columns, axis=1).sort_index().interpolate(limit_direction="both")
    bins = matrix.index.to_numpy()
    if len(bins) < 4:
        return bins * PENALTY_BIN, []
    result = getattr(rpt, algo_name)(model=model, min_size=2).fit(matrix.to_numpy()).predict(
        pen=pen * matrix.shape[1] * np.log(len(bins))
    )

    breakpoints = []
    for line_num, i in enumerate(result[:-1], start=1):
        lo = bins[max(0, i - REFINE_BINS)] * PENALTY_BIN
        hi = (bins[min(len(bins) - 1, i + REFINE_BINS)] + 1) * PENALTY_BIN
        refined = []
        for _bam_name, d, _span in sample_list:
            pos = d["base_pos" if "base_pos" in d.columns else "pos"].to_numpy()
            rows = np.flatnonzero((pos >= lo) & (pos < hi))
            if len(rows) >= 100:
                row = refine_breakpoint(d["norm"].to_numpy()[rows], 0, len(rows))
                refined.append(pos[rows[row]])
        breakpoints.append({
            "position": int(np.median(refined)) if refined else int(bins[i] * PENALTY_BIN),
            "breakpoint_number": line_num,
            "method": "ruptures_pooled_multivariate",
        })
    return bins * PENALTY_BIN, breakpoints


def run_ruptures_multivariate(
    sample_list: list[tuple[str, pd.DataFrame, str]],
    signal_col: str = "moving_average",
    algo_name: str = "BottomUp",
    model: str = "l2",
    n_bkps: int = 2,
    pen: float | None = None,
) -> tuple[np.ndarray, list[dict]]:
    """Run multivariate ruptures on a pooled set of samples.

    With *n_bkps*, samples are aligned to a shared position array (union of
    all positions) and interpolated to fill gaps before stacking into a
    signal matrix. With *pen*, see :func:`_pooled_penalty_breakpoints`.

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
    if pen is not None:
        return _pooled_penalty_breakpoints(sample_list, algo_name, model, pen)

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
