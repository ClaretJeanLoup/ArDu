"""
Signal-processing utilities for depth-of-coverage tracks.

Covers:
    - sliding-window median
    - sliding-window variance
    - sliding-window covariance (mean × variance, normalised to DoC scale)
    - rolling-average shift detection
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Sliding-window statistics
# ---------------------------------------------------------------------------

def window_median(arr: np.ndarray, w: int) -> pd.Series:
    """Sliding-window median over *arr* with window size *w*.

    Uses centred windows with ``min_periods=1`` so edge values are
    always defined.

    Parameters
    ----------
    arr:
        1-D array of coverage values (or any numeric signal).
    w:
        Window width in samples (not base-pairs).

    Returns
    -------
    pd.Series
        Same length as *arr*.
    """
    return pd.Series(arr).rolling(w, center=True, min_periods=1).median()


def window_variance(arr: np.ndarray, w: int) -> pd.Series:
    """Sliding-window variance over *arr* with window size *w*.

    Parameters
    ----------
    arr:
        1-D array of coverage values.
    w:
        Window width in samples.

    Returns
    -------
    pd.Series
        Same length as *arr*.
    """
    return pd.Series(arr).rolling(w, center=True, min_periods=1).var()


# ---------------------------------------------------------------------------
# Covariance track
# ---------------------------------------------------------------------------

def compute_sliding_covariance(
    data: pd.DataFrame,
    window_size: int,
) -> tuple[pd.Series, pd.Series]:
    """Compute covariance of sliding mean depth and sliding variance.

    The result is re-scaled to the same range as the ``"norm"`` column and
    shifted so that its mean equals 1; values below 1 are set to NaN.

    Parameters
    ----------
    data:
        DataFrame that must contain columns ``"depth"``, ``"norm"``, and
        ``"pos"``.
    window_size:
        Rolling-window width in samples.

    Returns
    -------
    pos_series : pd.Series
        Positional column (passthrough of ``data["pos"]``).
    cov_series : pd.Series
        Re-scaled covariance values aligned to the DoC axis.
    """
    mean_series = data["depth"].rolling(window_size, center=True).mean()
    var_series = data["depth"].rolling(window_size, center=True).var()
    cov_series = mean_series.rolling(window_size, center=True).cov(var_series)

    if not cov_series.isna().all():
        cov_min, cov_max = cov_series.min(), cov_series.max()
        norm_min, norm_max = data["norm"].min(), data["norm"].max()
        if cov_max != cov_min:  # avoid division by zero
            cov_series = (
                (cov_series - cov_min) / (cov_max - cov_min)
                * (norm_max - norm_min)
                + norm_min
            )
        cov_series = cov_series - cov_series.mean() + 1
        cov_series[cov_series < 1] = float("nan")

    return data["pos"], cov_series


# ---------------------------------------------------------------------------
# Shift detection
# ---------------------------------------------------------------------------

def detect_shifts(
    data: pd.DataFrame,
    threshold: float,
    window_size: int,
    passes: int,
) -> pd.DataFrame:
    """Detect positions where the moving average departs significantly from
    the normalised coverage.

    Both signals are smoothed *passes* times with a median filter of
    *window_size* before computing their difference.  Rows where the
    absolute difference exceeds *threshold* are returned.

    Parameters
    ----------
    data:
        DataFrame with columns ``"norm"`` and ``"moving_average"``.
    threshold:
        Absolute difference above which a position is flagged.
    window_size:
        Sliding-window width for additional smoothing.
    passes:
        Number of smoothing iterations.

    Returns
    -------
    pd.DataFrame
        Subset of *data* at flagged positions.
    """
    norm_smooth = data["norm"].values.copy()
    ma_smooth = data["moving_average"].values.copy()

    for _ in range(passes):
        norm_smooth = window_median(norm_smooth, window_size)
        ma_smooth = window_median(ma_smooth, window_size)

    diff = np.asarray(ma_smooth) - np.asarray(norm_smooth)
    return data[np.abs(diff) > threshold]


def bootstrap_median_ratio(
    num: np.ndarray,
    den: np.ndarray,
    n_boot: int,
    block: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    """95% interval of ``median(num) / median(den)`` by block bootstrap.

    Neighbouring bases share reads, so bases are resampled in contiguous
    blocks of *block* bp, independently in the locus and the normalisation
    locus. Returns ``(nan, nan)`` when either array is empty.
    """
    if num.size == 0 or den.size == 0:
        return float("nan"), float("nan")

    def blocks(a):
        b = max(1, min(block, a.size))
        return [a[i:i + b] for i in range(0, a.size, b)]

    nb, db = blocks(num), blocks(den)
    ratios = np.empty(n_boot)
    for k in range(n_boot):
        n = np.concatenate([nb[i] for i in rng.integers(0, len(nb), len(nb))])
        d = np.concatenate([db[i] for i in rng.integers(0, len(db), len(db))])
        md = np.median(d)
        ratios[k] = np.median(n) / md if md > 0 else np.nan
    ratios = ratios[~np.isnan(ratios)]
    if ratios.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(ratios, 2.5)), float(np.percentile(ratios, 97.5))
