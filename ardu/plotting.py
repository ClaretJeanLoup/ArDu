"""
Matplotlib plotting functions.

Covers:
    - per-sample depth-of-coverage plot
    - pooled multi-sample plot with optional grouped ruptures
    - parameter summary text box helper
    - plot-format output helpers
"""

from __future__ import annotations

import os
from collections import defaultdict

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from ardu.breakpoints import run_ruptures_multivariate
from ardu.probe import cluster_probe_profiles


matplotlib.use("Agg")  # non-interactive backend; safe for multiprocessing


# ---------------------------------------------------------------------------
# Parameter summary box
# ---------------------------------------------------------------------------

def add_plot_params_box(fig: matplotlib.figure.Figure, args) -> None:
    """Add a summary of run parameters inside the figure (lower-left).

    Parameters
    ----------
    fig:
        The target matplotlib Figure.
    args:
        Parsed argparse Namespace.
    """
    lines = []
    if getattr(args, "plot_threshold", None) is not None:
        lines.append(f"Threshold: {args.plot_threshold}")
    if getattr(args, "plot_slw", None) is not None:
        lines.append(f"Sliding Window: {args.plot_slw} bp")
    if getattr(args, "plot_bin", None) is not None:
        lines.append(f"Bin: {args.plot_bin} bp")
    if getattr(args, "plot_ylim", None):
        lines.append(f"Y-limits: {args.plot_ylim[0]}–{args.plot_ylim[1]}")
    if getattr(args, "plot_doclim", None):
        lines.append(f"DoC limits: {args.plot_doclim[0]}–{args.plot_doclim[1]}")
    if getattr(args, "plot_force", False):
        lines.append("Force Plotting: ON")
    fig.text(
        0.05,
        0.025,
        "\n".join(lines),
        fontsize=8,
        va="top",
        ha="left",
        bbox=dict(boxstyle="square", facecolor="white", edgecolor="gray", alpha=0.9),
    )


# ---------------------------------------------------------------------------
# Per-sample plot
# ---------------------------------------------------------------------------

def render_per_sample_plot(
    bam_name: str,
    locus: str,
    d: pd.DataFrame,
    str_incr_span: str,
    args,
    output_dir: str,
    breakpoints_data: list[dict],
    regions_dict: dict,
) -> None:
    """Render and save a depth-of-coverage plot for a single sample + locus.

    Parameters
    ----------
    bam_name:
        Sample identifier (BAM basename without extension).
    locus:
        Locus name.
    d:
        Per-base depth DataFrame with columns ``pos``, ``norm``,
        ``moving_average``, ``moving_variances``, ``covariance``.
    str_incr_span:
        SAM-style region string used as x-axis label.
    args:
        Parsed argparse Namespace.
    output_dir:
        Root output directory; a subdirectory ``{locus}-plots/`` is
        created inside it.
    breakpoints_data:
        List to append breakpoint dicts to (mutated in-place).
    regions_dict:
        Pre-parsed regions dict (for ``--plot-target`` highlight).
    """
    from ardu.breakpoints import run_ruptures_univariate, run_rollingaverage
    from ardu.intervals import get_totalspan

    fig, ax = plt.subplots(figsize=(11, 7))
    if args.plot_ylim:
        ax.set_ylim(args.plot_ylim)
    ax.grid(axis="y", linestyle="dashdot", linewidth=0.5)
    ax.tick_params(axis="both", labelsize=12)
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    ax.plot(
        d.pos, d.norm,
        label="Normalised Depth of Coverage",
        color="lightsteelblue",
        linestyle="none",
        marker=".",
    )

    if getattr(args, "plot_covar", False):
        ax.plot(d.pos, d.covariance, label="Sliding Window Covariance",
                color="blue", linewidth=2)
        ax.plot(d.pos, d.moving_variances, label="Sliding Window Variance",
                color="red", alpha=0.5, linewidth=2)
        ax.plot(d.pos, d.moving_average, label="Smoothed Depth of Coverage",
                color="black", alpha=0.5, linewidth=2)
    else:
        ax.plot(d.pos, d.moving_average, label="Smoothed Depth of Coverage",
                color="black", linewidth=2)

    # ------------------------------------------------------------------
    # Breakpoints
    # ------------------------------------------------------------------
    result = None
    regions_bkp = []

    if getattr(args, "breakpoint", None) == "ruptures":
        signal_col = args.bkp_signal if args.bkp_signal else "moving_average"
        bkps = run_ruptures_univariate(
            d,
            algo_name=args.bkp_algo,
            model=args.bkp_model,
            signal_col=signal_col,
            n_bkps=args.bkp_nb,
            pen=args.bkp_pen,
        )
        result = bkps

        x_range = d.pos.max() - d.pos.min()
        y_positions = []
        base_y = 1.4
        jitter_step = 0.15

        for entry in bkps:
            b = entry["row_index"]
            line_num = entry["breakpoint_number"]
            try:
                bp_pos = int(d.loc[b, "pos"])
            except (KeyError, ValueError, TypeError):
                print(f"Warning: invalid position at row {b}, skipping.")
                continue
            breakpoints_data.append(
                {
                    "bam_file": bam_name,
                    "locus": locus,
                    "chromosome": str_incr_span.split(":")[0],
                    "position": bp_pos,
                    "breakpoint_number": line_num,
                    "method": "ruptures",
                }
            )
            ax.axvline(bp_pos, color="darkred", linestyle="-", linewidth=1.5)

            y = base_y
            for prev_pos, prev_y in y_positions:
                if abs(bp_pos - prev_pos) < x_range * 0.05:
                    y = prev_y + jitter_step
            y_positions.append((bp_pos, y))

            ax.text(
                bp_pos, y, f"{line_num}:{bp_pos}",
                color="darkred", fontsize=10, weight="bold", va="bottom",
            )

    elif getattr(args, "breakpoint", None) == "rollingaverage":
        regions_bkp = run_rollingaverage(
            d,
            threshold=args.bkp_threshold,
            window_size=args.bkp_slw,
            passes=args.bkp_passes,
        )
        for i, region in enumerate(regions_bkp):
            region_center = (region[0] + region[1]) / 2
            region_data = d[(d["pos"] >= region[0]) & (d["pos"] <= region[1])]
            if not region_data.empty:
                max_norm_depth = region_data.loc[region_data["norm"].idxmax(), "norm"]
                ax.axvspan(region[0], region[1], color="darkred", alpha=0.3)
                ax.text(
                    region_center, max_norm_depth + 0.1, str(i + 1),
                    color="darkred", fontsize=10, weight="bold",
                )

    # ------------------------------------------------------------------
    # Optional target highlight
    # ------------------------------------------------------------------
    if getattr(args, "plot_target", False):
        totalspan = get_totalspan(regions_dict, locus)
        if totalspan:
            mask = (d.pos >= totalspan[1]) & (d.pos <= totalspan[2])
            ax.plot(
                d.pos[mask], d.moving_average[mask],
                color="mediumpurple", label=locus, zorder=5, linewidth=2.5,
            )

    # ------------------------------------------------------------------
    # Labels, legend, save
    # ------------------------------------------------------------------
    fig.suptitle(f"{locus} locus in {bam_name}", fontsize=14, y=1.02)
    ax.set_xlabel(f"Genomic position: {str_incr_span}", fontsize=12)
    ax.set_ylabel("Normalised Depth of Coverage", fontsize=12)

    if getattr(args, "breakpoint", None):
        custom_line = Line2D(
            [0], [0], color="darkred", linestyle="-", linewidth=1.5,
            label="Predicted Breakpoints",
        )
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles + [custom_line], labels + ["Predicted Breakpoints"],
                  loc="upper right")
    else:
        ax.legend(loc="upper right")

    if getattr(args, "plot_param", False):
        plt.subplots_adjust(bottom=0.25)
        add_plot_params_box(fig, args)

    plt.tight_layout()
    output_dirlocus = os.path.join(output_dir, f"{locus}-plots/")
    os.makedirs(output_dirlocus, exist_ok=True)
    plt.savefig(
        os.path.join(output_dirlocus, f"{bam_name}_{locus}_plot.{args.plot}"),
        bbox_inches="tight",
        dpi=args.plot_dpi,
    )
    plt.close()


# ---------------------------------------------------------------------------
# Pooled plot
# ---------------------------------------------------------------------------

def make_pooled_plot(
    locus: str,
    sample_list: list[tuple[str, pd.DataFrame, str]],
    args,
    output_dir: str,
    breakpoints_data: list[dict] | None = None,
    probe_profiles: dict[str, np.ndarray] | None = None,
    group_label: str = "pooled",
) -> None:
    """Render a pooled depth-of-coverage plot for *locus*.

    When *probe_profiles* is supplied, samples are first grouped by
    structural similarity via :func:`ardu.probe.cluster_probe_profiles`
    and one plot is produced per group.  When ``--breakpoint ruptures`` is
    active, multivariate ruptures runs once per group on the joint signal
    matrix so that breakpoints reflect positions consistent across all
    samples in that group.

    Parameters
    ----------
    locus:
        Locus name (used in figure title and filename).
    sample_list:
        List of ``(bam_name, depth_df, span_str)`` tuples.
    args:
        Parsed argparse Namespace.
    output_dir:
        Root output directory.
    breakpoints_data:
        List to append breakpoint dicts to (mutated in-place).  Pass
        ``None`` to discard breakpoint records.
    probe_profiles:
        ``{bam_name: ratio_array}`` — when present, triggers per-group
        plotting.  Pass ``None`` to plot all samples together.
    group_label:
        String suffix appended to the output filename
        (e.g. ``"pooled"``, ``"group0"``, ``"group1"``).
    """
    # ------------------------------------------------------------------
    # Probe-profile grouping: recurse per group
    # ------------------------------------------------------------------
    if probe_profiles and len(probe_profiles) > 1:
        labels = cluster_probe_profiles(
            probe_profiles,
            threshold=getattr(args, "probe_threshold", 1.4),
            distance_threshold=getattr(args, "probe_cluster_dist", 0.15),
        )
        groups: dict[int, list] = defaultdict(list)
        for entry in sample_list:
            bam_name = entry[0]
            groups[labels.get(bam_name, 0)].append(entry)

        print(f"  {locus}: {len(groups)} probe-profile group(s) identified")
        for gid, group_members in sorted(groups.items()):
            member_names = [m[0] for m in group_members]
            print(f"    Group {gid}: {member_names}")
            make_pooled_plot(
                locus,
                group_members,
                args,
                output_dir,
                breakpoints_data,
                probe_profiles=None,  # don't re-cluster within a group
                group_label=f"group{gid}",
            )
        return

    # ------------------------------------------------------------------
    # Single pooled plot
    # ------------------------------------------------------------------
    cmap = plt.get_cmap("tab10")
    fig, ax = plt.subplots(figsize=(13, 7))
    if getattr(args, "plot_ylim", None):
        ax.set_ylim(args.plot_ylim)
    ax.grid(axis="y", linestyle="dashdot", linewidth=0.5)
    ax.tick_params(axis="both", labelsize=12)
    plt.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

    for i, (bam_name, d, _span) in enumerate(sample_list):
        color = cmap(i % 10)
        if getattr(args, "plot_pooled_raw", False):
            ax.plot(d.pos, d.norm, color=color, linestyle="none",
                    marker=".", alpha=0.1, markersize=2)
        ax.plot(d.pos, d.moving_average, label=bam_name, color=color,
                linewidth=1.5, alpha=0.85)

    # ------------------------------------------------------------------
    # Multivariate ruptures
    # ------------------------------------------------------------------
    bkp_result = None
    if getattr(args, "breakpoint", None) == "ruptures" and sample_list:
        signal_col = args.bkp_signal if args.bkp_signal else "moving_average"
        span_str = sample_list[0][2]
        ref_d = sample_list[0][1]
        x_range = ref_d.pos.max() - ref_d.pos.min()

        _ref_pos, bkp_list = run_ruptures_multivariate(
            sample_list,
            signal_col=signal_col,
            algo_name=args.bkp_algo,
            model=args.bkp_model,
            n_bkps=args.bkp_nb,
            pen=args.bkp_pen,
        )
        bkp_result = bkp_list

        y_positions = []
        base_y = 1.4
        jitter_step = 0.15

        for entry in bkp_list:
            bp_pos = entry["position"]
            line_num = entry["breakpoint_number"]
            ax.axvline(bp_pos, color="darkred", linestyle="--", linewidth=1.5)

            y = base_y
            for prev_pos, prev_y in y_positions:
                if abs(bp_pos - prev_pos) < x_range * 0.05:
                    y = prev_y + jitter_step
            y_positions.append((bp_pos, y))

            ax.text(
                bp_pos, y, f"{line_num}:{bp_pos}",
                color="darkred", fontsize=10, weight="bold", va="bottom",
            )
            if breakpoints_data is not None:
                breakpoints_data.append(
                    {
                        "bam_file": "pooled",
                        "locus": locus,
                        "chromosome": span_str.split(":")[0],
                        "position": bp_pos,
                        "breakpoint_number": line_num,
                        "method": "ruptures_pooled_multivariate",
                    }
                )

    # ------------------------------------------------------------------
    # Legend, title, save
    # ------------------------------------------------------------------
    if bkp_result:
        custom_line = Line2D(
            [0], [0], color="darkred", linestyle="--", linewidth=1.5,
            label="Pooled Breakpoints (multivariate)",
        )
        handles, labels = ax.get_legend_handles_labels()
        ax.legend(
            handles + [custom_line],
            labels + ["Pooled Breakpoints (multivariate)"],
            loc="upper right",
            fontsize=7,
            ncol=max(1, len(sample_list) // 15),
        )
    else:
        ax.legend(
            loc="upper right",
            fontsize=7,
            ncol=max(1, len(sample_list) // 15),
        )

    fig.suptitle(
        f"{locus} — {group_label} (n={len(sample_list)})",
        fontsize=14, y=1.02,
    )
    ax.set_xlabel("Genomic position", fontsize=12)
    ax.set_ylabel("Normalised Depth of Coverage", fontsize=12)

    if getattr(args, "plot_param", False):
        plt.subplots_adjust(bottom=0.25)
        add_plot_params_box(fig, args)

    plt.tight_layout()
    output_dirlocus = os.path.join(output_dir, f"{locus}-plots/")
    os.makedirs(output_dirlocus, exist_ok=True)
    suffix = f"_{group_label}_ruptures" if bkp_result else f"_{group_label}"
    plt.savefig(
        os.path.join(output_dirlocus, f"{locus}{suffix}.{args.plot}"),
        bbox_inches="tight",
        dpi=args.plot_dpi,
    )
    plt.close()