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
