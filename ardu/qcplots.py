"""Diagnostic figures of ``ardu coverage``: zero-depth bases per locus, and GC against RCN spread."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

BLUE, BLUE_MID, BLUE_LIGHT = "#2a78d6", "#7aa9e8", "#c9ddf6"
ORANGE, GREY, INK, MUTED = "#eb6834", "#a3a29d", "#0b0b0b", "#52514e"
STYLE = {
    "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": MUTED, "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK,
    "axes.titlesize": 10, "axes.titleweight": "bold", "axes.titlelocation": "left",
}


def plot_zero_depth(region_df: pd.DataFrame, path: str) -> pd.DataFrame:
    """Per locus: % of bases at depth 0 in each sample, and samples above 0, 10 and 50%.

    *region_df* has the columns of ``_region_qc.tsv`` (``bam_file``, ``locus``, ``length``,
    ``pct_zero``). Returns the per-locus summary drawn, sorted as on the figure.
    """
    df = region_df[pd.to_numeric(region_df["pct_zero"], errors="coerce").notna()].copy()
    df["zero_bases"] = df["length"] * df["pct_zero"].astype(float) / 100
    per = df.groupby(["locus", "bam_file"], sort=False)[["length", "zero_bases"]].sum()
    per["pct"] = 100 * per["zero_bases"] / per["length"]
    pct = per["pct"].unstack()                       # loci x samples
    summary = pd.DataFrame({
        "length": per["length"].groupby(level=0).first(),
        "median": pct.median(axis=1),
        "mean": pct.mean(axis=1),
        "gt0": (pct > 0).sum(axis=1),
        "ge10": (pct >= 10).sum(axis=1),
        "ge50": (pct >= 50).sum(axis=1),
    }).sort_values(["median", "mean"], ascending=False)
    loci = list(summary.index)
    n_samples = pct.shape[1]
    # zoom on the data, so that loci all under a few % stay readable; never below 10% or above 100%
    xmax = min(102.0, max(10.0, float(np.nanmax(pct.to_numpy())) * 1.1))

    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(1, 2, figsize=(13, 0.28 * len(loci) + 1.6), sharey=True,
                               gridspec_kw={"width_ratios": [1.6, 1]})
        a = ax[0]
        for i, locus in enumerate(loci):
            v = pct.loc[locus].dropna().to_numpy()
            # golden-ratio jitter: spread without randomness, so reruns give the same figure
            jitter = i + ((np.arange(v.size) * 0.6180339) % 1 - 0.5) * 0.6
            a.scatter(v, jitter, s=6, c=GREY, alpha=0.6, linewidths=0)
            a.plot([np.median(v)], [i], marker="|", ms=14, mew=2.5, color=BLUE)
        a.set(yticks=range(len(loci)),
              yticklabels=[f"{locus} ({int(summary.loc[locus, 'length'])} bp)" for locus in loci],
              xlim=(-xmax / 50, xmax), xlabel="% of locus bases at depth 0 (one dot per sample, bar = median)",
              title="A  zero-depth bases per locus")
        a.invert_yaxis()
        a.grid(axis="x", color="#e6e6e3", lw=0.8)
        a.set_axisbelow(True)

        a = ax[1]
        h = 0.27
        for j, (col, label, colour) in enumerate((("gt0", "> 0%", BLUE_LIGHT), ("ge10", ">= 10%", BLUE_MID),
                                                  ("ge50", ">= 50%", BLUE))):
            a.barh([i + (j - 1) * h for i in range(len(loci))], summary[col], height=h - 0.03,
                   color=colour, label=f"{label} at depth 0")
        a.set(xlabel=f"samples (of {n_samples})", xlim=(0, n_samples * 1.02),
              title="B  samples with zero-depth bases")
        a.xaxis.set_major_locator(MaxNLocator(integer=True))
        a.legend(frameon=False, loc="lower right")
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return summary


def plot_gc_rcn(rcn: pd.DataFrame, gc_by_locus: dict, norm: str, gc_warn: float, path: str) -> pd.DataFrame:
    """Per locus: GC difference from the normalisation against the spread of RCN across samples.

    *rcn* is loci x samples. The spread is the SD of log2 RCN over the samples with RCN > 0;
    a locus far from the normalisation in GC tends to have a noisier RCN, as does a real CNV.
    """
    norm_gc = gc_by_locus[norm][0]
    rows = []
    for locus, values in rcn.iterrows():
        if locus == norm or locus not in gc_by_locus or np.isnan(gc_by_locus[locus][0]):
            continue
        v = pd.to_numeric(values, errors="coerce")
        v = v[v > 0]
        if v.size < 2:
            continue
        rows.append((locus, 100 * (gc_by_locus[locus][0] - norm_gc), float(np.std(np.log2(v), ddof=1)), v.size))
    table = pd.DataFrame(rows, columns=["locus", "gc_diff", "log2_sd", "n"]).set_index("locus")
    if table.empty:
        return table

    far = table["gc_diff"].abs() > gc_warn
    with plt.rc_context(STYLE):
        fig, a = plt.subplots(figsize=(8, 5.5))
        a.axvspan(-gc_warn, gc_warn, color="#e6efec", zorder=0)
        a.axvline(0, color="#7a7f88", lw=1)
        a.scatter(table.loc[~far, "gc_diff"], table.loc[~far, "log2_sd"], s=60, c="#33373d", alpha=0.75,
                  edgecolors="white", linewidths=2, zorder=3)
        a.scatter(table.loc[far, "gc_diff"], table.loc[far, "log2_sd"], s=60, c=ORANGE,
                  edgecolors="white", linewidths=2, zorder=3, label=f"more than {gc_warn:g} GC points away")
        # every locus is named when they are few, otherwise the GC-distant ones and the noisiest
        named = table.index if len(table) <= 25 else table.index[far | (table["log2_sd"] >= table["log2_sd"].quantile(0.9))]
        for locus in named:
            a.annotate(locus, (table.loc[locus, "gc_diff"], table.loc[locus, "log2_sd"]), xytext=(6, 4),
                       textcoords="offset points", fontsize=8, color=INK)
        lim = max(table["gc_diff"].abs().max(), gc_warn) * 1.15
        top = table["log2_sd"].max() * 1.18
        a.text(0, top * 0.98, f"within {gc_warn:g} GC points", ha="center", va="top", color=MUTED)
        a.set(xlim=(-lim, lim), ylim=(0, top), xlabel=f"GC difference from {norm} (percentage points)",
              ylabel="spread of RCN across samples (SD of log2 RCN)",
              title=f"GC distance to the normalisation and RCN noise ({int(table['n'].max())} samples)")
        a.grid(axis="y", color="#e6e6e3", lw=0.8)
        a.set_axisbelow(True)
        if far.any():
            a.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.12))
        fig.tight_layout()
        fig.savefig(path, dpi=200)
        plt.close(fig)
    return table
