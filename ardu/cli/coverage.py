"""
``ardu coverage`` — depth-of-coverage analysis and breakpoint detection.

Corresponds to the original ArDu_1_1.py, refactored into:
    - register(subparsers)  – adds the 'coverage' subcommand
    - run(args)             – orchestration loop
    - helper functions      – _resolve_plot_interval, _build_depth_df,
                              _process_one_bam, _run_parallel
"""

from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import pysam

from ardu.bam import (
    calculate_coverage_stats,
    calculate_indel_counts,
    configure_read_filters,
    calculate_nucleotide_counts,
    depth_over_region,
    parse_region,
)
from ardu.gc import (
    GC_PAD,
    draw_norm_windows,
    loci_gc,
    read_bed_intervals,
    write_gc_report,
    write_windows_bed,
)
from ardu.intervals import get_regions, get_totalspan, parse_plot_intervals
from ardu.plotting import make_pooled_plot, render_per_sample_plot
from ardu.probe import build_probe_profile, write_probe_profiles_tsv
from ardu.signal import bootstrap_median_ratio, compute_sliding_covariance, window_median, window_variance


# ---------------------------------------------------------------------------
# Argparse registration
# ---------------------------------------------------------------------------

def register(subparsers) -> None:
    """Add the ``coverage`` sub-command to *subparsers*."""
    parser = subparsers.add_parser(
        "coverage",
        help="Depth-of-coverage analysis, normalisation, and breakpoint detection.",
        description=(
            "ArDu coverage — identify duplication structure via depth of "
            "coverage analysis."
        ),
    )

    # Mandatory
    parser.add_argument("-b", "--bam", required=True,
                        help="File listing BAM paths, one per line.")
    parser.add_argument("-r", "--region", required=True,
                        help="Four-column BED-like file "
                             "(chromosome, start, stop, locus_name).  "
                             "Regions sharing a locus_name are pooled.")
    norm = parser.add_mutually_exclusive_group(required=True)
    norm.add_argument("-n", "--norm",
                      help="Name of the normalisation locus in the region file, or a region file "
                           "of its own (e.g. the _norm_windows.bed of an earlier run).")
    norm.add_argument("--norm-random", type=int, metavar="N",
                      help="Normalise on N random windows instead of a locus (needs --fasta). "
                           "They are written to <prefix>_norm_windows.bed for reuse with -n.")
    parser.add_argument("-o", "--outfile", required=True,
                        help="Prefix for output file names. The run folder is created in "
                             "the prefix's directory (current directory if none).")

    # GC and random normalisation windows
    parser.add_argument("--fasta", metavar="REFERENCE",
                        help="Reference fasta (indexed) of the bams: writes <prefix>_gc.tsv and "
                             "warns when a locus and the normalisation differ in GC.")
    parser.add_argument("--gc-warn", type=float, default=3, metavar="POINTS",
                        help="GC difference (percentage points) above which a locus is flagged.")
    parser.add_argument("--norm-size", type=int, default=1000,
                        help="Size (bp) of the random normalisation windows.")
    parser.add_argument("--norm-gc", choices=["match", "any"], default="match",
                        help="'match': keep windows within 2 GC points of the median GC of the "
                             "loci; 'any': no GC selection.")
    parser.add_argument("--norm-ignore", nargs="+", metavar="LOCUS", default=[],
                        help="Loci of the region file left out of the GC the windows match "
                             "(e.g. an old reference locus).")
    parser.add_argument("--norm-exclude", metavar="BED",
                        help="BED of regions the windows must avoid (known CNVs, repeats).")
    parser.add_argument("--norm-chroms", nargs="+", metavar="CONTIG",
                        help="Contigs to draw windows from (e.g. the autosomes); "
                             "default: every contig of at least 1 Mb.")
    parser.add_argument("--norm-margin", type=int, default=50_000,
                        help="Distance (bp) kept between the windows and the loci of the region file.")
    parser.add_argument("--norm-seed", type=int, default=1,
                        help="Random seed of the window draw.")

    # Read filters
    parser.add_argument("--min-mapq", type=int, default=None,
                        help="Minimum mapping quality of reads counted in depth and "
                             "mutation counts. Default=20. Lower it for loci with "
                             "identical paralog copies in the reference, whose reads "
                             "get MAPQ 0.")
    parser.add_argument("--min-baseq", type=int, default=None,
                        help="Minimum base quality counted in depth and mutation "
                             "counts. Default=13.")

    # RCN uncertainty
    parser.add_argument("--rcn-ci", type=int, default=0, metavar="N",
                        help="Add a 95%% confidence interval to each RCN from N block "
                             "bootstrap replicates (two extra fields per cell). "
                             "Default=0 (off).")
    parser.add_argument("--ci-block", type=int, default=500,
                        help="Block size (bp) of the RCN bootstrap. Default=500.")
    parser.add_argument("--ci-seed", type=int, default=1,
                        help="Random seed of the RCN bootstrap. Default=1.")

    # Parallelism
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of parallel BAM processing workers. Default=1.")

    # Plotting
    parser.add_argument("--plot", type=str, default=None,
                        help="Produce per-sample plots. Accepts: png, jpeg, jpg, pdf, svg, eps.")
    parser.add_argument("--plot-threshold", type=float, default=1.4,
                        help="Minimum normalised coverage to trigger plotting. Default=1.4.")
    parser.add_argument("--plot-interval",
                        help="Tab-delimited file of custom plot intervals: locus<TAB>chrom:start-stop")
    parser.add_argument("--plot-proportion", type=float, default=2,
                        help="Extend plot region by this multiple of locus span. Default=2.")
    parser.add_argument("--plot-auto", action="store_true",
                        help="Automatically expand the plot interval by probing coverage.")
    parser.add_argument("--probe-size", type=int, default=500,
                        help="Probe window size (bp). Default=500.")
    parser.add_argument("--max-extension", type=int, default=5_000_000,
                        help="Maximum interval extension on each side (bp). Default=5 Mb.")
    parser.add_argument("--probe-threshold", type=float, default=0.8,
                        help="Coverage ratio of probes to target. Default=0.8.")
    parser.add_argument("--probe-number", type=int, default=20,
                        help="Number of probes per round of extension. Default=20.")
    parser.add_argument("--probe-spacing", type=int, default=1000,
                        help="Distance between probe starts (bp). Default=1000.")
    parser.add_argument("--probe-drops", type=int, default=10,
                        help="Stop after this many consecutive low-coverage probes. Default=10.")
    parser.add_argument("--probe-symmetric", action="store_true",
                        help="Enable symmetric extension (default: OFF).")
    parser.add_argument("--plot-slw", type=int, default=1000,
                        help="Sliding window size (bp) for coverage smoothing. Default=1000.")
    parser.add_argument("--plot-target", action="store_true",
                        help="Highlight the target position on the plot.")
    parser.add_argument("--plot-force", action="store_true",
                        help="Force plotting regardless of coverage value.")
    parser.add_argument("--plot-covar", action="store_true",
                        help="Overlay covariance and variance tracks.")
    parser.add_argument("--plot-ylim", type=float, nargs=2,
                        help="Y-axis limits: two numbers, e.g. 0 4.")
    parser.add_argument("--plot-doclim", type=float, nargs=2,
                        help="Min and max normalised DoC for plotting.")
    parser.add_argument("--plot-bin", type=int,
                        help="Genomic bin size (bp) applied on top of sliding window.")
    parser.add_argument("--plot-dpi", type=int, default=300,
                        help="Plot resolution (DPI). Default=300.")
    parser.add_argument("--plot-param", action="store_true",
                        help="Print run parameters on the plot.")

    # Pooled plotting
    parser.add_argument("--plot-pooled", action="store_true",
                        help="Produce an additional pooled plot per locus.  "
                             "Combined with --breakpoint ruptures, runs multivariate "
                             "ruptures on all samples jointly.")
    parser.add_argument("--plot-pooled-raw", action="store_true",
                        help="Show per-base normalised scatter on pooled plots.")

    # Breakpoints
    parser.add_argument("-bkp", "--breakpoint",
                        choices=["ruptures", "rollingaverage"],
                        help="Breakpoint detection method.")
    parser.add_argument("--bkp-slw", type=int, default=1000,
                        help="Window size for rolling-average breakpoint detection. Default=1000.")
    parser.add_argument("--bkp-signal", type=str,
                        help="Signal column for ruptures: norm, moving_average, "
                             "moving_variances, or covariance. Default=moving_average.")
    parser.add_argument("--bkp-nb", type=int, default=2,
                        help="Expected number of breakpoints (ruptures). Default=2.")
    parser.add_argument("--bkp-pen", type=float,
                        help="Find the number of breakpoints with a penalty instead of "
                             "--bkp-nb: ruptures runs on the mean normalised depth per bin "
                             "(--plot-bin, else 500 bp) with a penalty of this value times "
                             "the noise variance times log(number of bins). About 10 suits "
                             "depth data; lower finds smaller steps.")
    parser.add_argument("--bkp-model", type=str, default="l2",
                        help="ruptures cost model. Default=l2.")
    parser.add_argument("--bkp-algo", type=str, default="BottomUp",
                        help="ruptures algorithm. Default=BottomUp.")
    parser.add_argument("--bkp-threshold", type=float, default=0.5,
                        help="Shift threshold for rollingaverage. Default=0.5.")
    parser.add_argument("--bkp-passes", type=int, default=1,
                        help="Smoothing passes for rollingaverage. Default=1.")

    # Probe profile
    parser.add_argument("--probe-profile", action="store_true",
                        help="Build per-sample probe profiles around each locus.  "
                             "Combined with --plot-pooled, groups samples by structural "
                             "similarity before running pooled ruptures.")
    parser.add_argument("--probe-cluster-dist", type=float, default=0.15,
                        help="Maximum Hamming distance for probe-profile clustering. Default=0.15.")
    parser.add_argument("--probe-n-each-side", type=int, default=20,
                        help="Number of probes on each side of the target. Default=20.")
    parser.add_argument("--probe-profile-out", action="store_true",
                        help="Write _probe_profiles.tsv with per-sample profile strings.")

    # Mutation genotyping
    parser.add_argument("--mutation",
                        help="Tab-delimited file: chromosome, position, [mutation_name].  "
                             "Returns per-nucleotide read counts.")
    parser.add_argument("--mutation-indels", action="store_true",
                        help="With --mutation, also count reads with an insertion or a "
                             "deletion right after each position (VCF-style anchor) and "
                             "the length of the homopolymer run that follows it in each "
                             "read, e.g. ins=8;del=0;runs=T8:12,T9:3,T10:5.")

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_bound(val, normalised_min: float, normalised_max: float) -> float:
    """Resolve a ``--plot-doclim`` bound value (float, 'min', or 'max')."""
    if val == "min":
        return normalised_min
    if val == "max":
        return normalised_max
    return float(val)


def _resolve_plot_interval(
    bam_file: str,
    locus: str,
    mean_coverage: float,
    ref_median: float,
    args,
    regions_dict: dict,
    plot_intervals_cache: dict | None,
) -> tuple[str | None, list[tuple[int, int]]]:
    """Return ``(span_str, depth_data)`` for the given locus.

    Three strategies are tried in order:
    1. Custom interval file (``--plot-interval``)
    2. Automatic expansion (``--plot-auto``)
    3. Proportional extension (``--plot-proportion``, the default)

    Returns ``(None, [])`` on any failure so the caller can skip the locus.
    """
    from ardu.probe import expand_interval

    if plot_intervals_cache is not None:
        str_incr_span = plot_intervals_cache.get(locus)
        if str_incr_span is None:
            print(f"Warning: {locus} not found in plot interval file. Skipping.")
            return None, []
        try:
            depth_data = depth_over_region(bam_file, str_incr_span)
        except (ValueError, KeyError) as exc:
            print(f"Warning: bad plot interval '{str_incr_span}' for {locus} ({exc}). Skipping.")
            return None, []
        return str_incr_span, depth_data

    if args.plot_auto:
        totalspan = get_totalspan(regions_dict, locus)
        if totalspan is None:
            return None, []
        chrom, dup_start, dup_end = totalspan
        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            chrom_lengths = dict(zip(samfile.references, samfile.lengths))
        chrom_len = chrom_lengths.get(chrom)
        if chrom_len is None:
            print(f"Warning: chromosome {chrom} not found in BAM header, skipping.")
            return None, []
        if dup_start >= dup_end:
            print(f"Warning: invalid interval {chrom}:{dup_start}-{dup_end}, skipping.")
            return None, []
        dup_start = max(0, dup_start)
        dup_end = min(chrom_len, dup_end)
        str_incr_span = expand_interval(
            bam_file, chrom, dup_start, dup_end,
            target_cov=mean_coverage,
            reference_cov=ref_median,
            probe_size=args.probe_size,
            probe_spacing=args.probe_spacing,
            probe_number=args.probe_number,
            stop_drops=args.probe_drops,
            max_extension=args.max_extension,
            symmetric=args.probe_symmetric,
        )
        depth_data = depth_over_region(bam_file, str_incr_span)
        return str_incr_span, depth_data

    if args.plot_proportion:
        totalspan = get_totalspan(regions_dict, locus)
        if totalspan is None:
            print(f"Warning: No intervals found for locus {locus}. Skipping.")
            return None, []
        chromosome, plotStart, plotStop = totalspan
        locus_length = plotStop - plotStart
        extension = locus_length * args.plot_proportion
        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            chrom_lengths = dict(zip(samfile.references, samfile.lengths))
        chrom_len = chrom_lengths.get(chromosome)
        try:
            plotStop = min(chrom_len, plotStop + int(extension))
            plotStart = max(1, plotStart - int(extension))
        except (ValueError, TypeError):
            print(f"Warning: invalid extension for locus {locus}. Skipping.")
            return None, []
        str_incr_span = f"{chromosome}:{plotStart}-{plotStop}"
        depth_data = depth_over_region(bam_file, str_incr_span)
        return str_incr_span, depth_data

    return None, []


def _build_depth_df(
    depth_data: list[tuple[int, int]],
    ref_median: float,
    args,
) -> pd.DataFrame:
    """Build the per-base depth DataFrame from raw pileup data.

    Applies optional DoC limits, computes sliding-window signals, and
    optionally bins the data.

    Parameters
    ----------
    depth_data:
        List of ``(position, depth)`` tuples from pysam pileup.
    ref_median:
        Normalisation constant.
    args:
        Parsed argparse Namespace.

    Returns
    -------
    pd.DataFrame
        Columns: ``pos``, ``depth``, ``norm``, ``moving_average``,
        ``moving_variances``, ``covariance``.
    """
    d = pd.DataFrame(depth_data, columns=["pos", "depth"])
    d["depth"] = pd.to_numeric(d["depth"], errors="coerce")
    d["norm"] = d["depth"] / ref_median

    if getattr(args, "plot_doclim", None):
        norm_min = d["norm"].min()
        norm_max = d["norm"].max()
        a = _parse_bound(args.plot_doclim[0], norm_min, norm_max)
        b = _parse_bound(args.plot_doclim[1], norm_min, norm_max)
        d = d[(d["norm"] >= a) & (d["norm"] <= b)].reset_index(drop=True)

    # binning below averages pos within bins; breakpoints are reported per base
    d["base_pos"] = d["pos"]
    d["moving_average"] = window_median(d["norm"].to_numpy(), args.plot_slw)
    d["moving_variances"] = window_variance(d["norm"].to_numpy(), args.plot_slw)
    d["pos"], d["covariance"] = compute_sliding_covariance(d, args.plot_slw)

    bin_size = getattr(args, "plot_bin", None)
    if bin_size and not d.empty:
        region_span = d["pos"].max() - d["pos"].min()
        if bin_size < region_span:
            bin_edges = np.arange(d["pos"].min(), d["pos"].max() + bin_size, bin_size)
            if len(bin_edges) > 1:
                d["genomic_bins"] = pd.cut(d["pos"], bins=bin_edges, include_lowest=True)
                for col in ["pos", "moving_average", "moving_variances", "covariance"]:
                    d[col] = d.groupby("genomic_bins", observed=True)[col].transform("mean")
            else:
                print(f"Skipping binning: not enough bins for bin_size={bin_size}.")
        else:
            print(f"Skipping binning: bin size {bin_size} >= region span {region_span}.")

    return d


# ---------------------------------------------------------------------------
# Per-BAM worker
# ---------------------------------------------------------------------------

def _process_one_bam(
    args_tuple: tuple,
) -> tuple[str, pd.DataFrame | None, list[dict], dict]:
    """Process a single BAM file.

    Designed to be called either directly in the serial loop or submitted
    to a ``ProcessPoolExecutor``.  Must be a top-level function for pickle
    compatibility.

    Parameters
    ----------
    args_tuple:
        ``(bam_file, args, regions_dict, plot_intervals_cache, output_dir)``

    Returns
    -------
    bam_name : str
    cov_df : pd.DataFrame | None
    breakpoints_list : list[dict]
    pooled_depth_dict : dict
        ``{locus: [(bam_name, depth_df, span_str), ...]}``
    probe_profiles_dict : dict
        ``{locus: {bam_name: ratio_array}}``
    """
    bam_file, args, regions_dict, plot_intervals_cache, output_dir = args_tuple
    configure_read_filters(args.min_mapq, args.min_baseq)
    bam_name = os.path.basename(bam_file).removesuffix(".bam")
    print(f"Processing {bam_name}")

    depth_arrays = {} if args.rcn_ci > 0 else None
    locus_coverage_df = calculate_coverage_stats(regions_dict, bam_file, depth_arrays)

    if args.norm not in locus_coverage_df["locus"].values:
        print(f"Normalisation locus not found in {bam_file}. Skipping.")
        return bam_name, None, [], {}, {}

    ref_median = locus_coverage_df.loc[
        locus_coverage_df["locus"] == args.norm, "median"
    ].values[0]
    if not isinstance(ref_median, (int, float)):
        print(f"Invalid normalisation coverage in {bam_file}. Skipping.")
        return bam_name, None, [], {}, {}

    locus_coverage_df["normalised"] = locus_coverage_df["median"].apply(
        lambda x: x / ref_median if isinstance(x, (int, float)) else "NA"
    )
    if depth_arrays is not None:
        rng = np.random.default_rng(args.ci_seed)
        bounds = [
            bootstrap_median_ratio(depth_arrays.get(locus, np.zeros(0)), depth_arrays[args.norm],
                                   args.rcn_ci, args.ci_block, rng)
            for locus in locus_coverage_df["locus"]
        ]
        locus_coverage_df["rcn_low"] = [b[0] for b in bounds]
        locus_coverage_df["rcn_high"] = [b[1] for b in bounds]
    locus_coverage_df["bam_file"] = bam_name

    if not args.plot:
        return bam_name, locus_coverage_df, [], {}, {}

    # Coerce numeric columns; warn about and drop NA rows
    for col in ["mean", "sd", "normalised"]:
        locus_coverage_df[col] = pd.to_numeric(locus_coverage_df[col], errors="coerce")

    na_rows = locus_coverage_df[
        locus_coverage_df[["mean", "sd", "normalised"]].isna().any(axis=1)
    ]
    if not na_rows.empty:
        print("Warning: skipping rows with missing values:")
        for locus_name in na_rows["locus"]:
            print(f"  - {locus_name}")
        locus_coverage_df = locus_coverage_df.dropna(
            subset=["mean", "sd", "normalised"]
        ).reset_index(drop=True)

    if args.plot_force:
        filtered_df = locus_coverage_df[locus_coverage_df["locus"] != args.norm]
    else:
        filtered_df = locus_coverage_df[
            locus_coverage_df["normalised"] > args.plot_threshold
        ]

    breakpoints_list: list[dict] = []
    pooled_depth_dict: dict = {}
    probe_profiles_dict: dict = {}

    for _idx, row in filtered_df.iterrows():
        locus = row["locus"]
        mean_coverage = row["mean"]

        str_incr_span, depth_data = _resolve_plot_interval(
            bam_file, locus, mean_coverage, ref_median,
            args, regions_dict, plot_intervals_cache,
        )
        if str_incr_span is None:
            continue

        d = _build_depth_df(depth_data, ref_median, args)

        # Probe profile
        if getattr(args, "probe_profile", False):
            totalspan = get_totalspan(regions_dict, locus)
            if totalspan:
                chrom_p, start_p, end_p = totalspan
                _, probe_ratios = build_probe_profile(
                    bam_file, chrom_p, start_p, end_p,
                    ref_median=ref_median,
                    probe_size=args.probe_size,
                    probe_spacing=args.probe_spacing,
                    n_probes_each_side=args.probe_n_each_side,
                )
                probe_profiles_dict.setdefault(locus, {})[bam_name] = probe_ratios

        # Store for pooled plot
        if getattr(args, "plot_pooled", False):
            pooled_depth_dict.setdefault(locus, []).append(
                (bam_name, d.copy(), str_incr_span)
            )

        # Per-sample plot
        render_per_sample_plot(
            bam_name, locus, d, str_incr_span, args, output_dir,
            breakpoints_list, regions_dict,
        )

    return bam_name, locus_coverage_df, breakpoints_list, pooled_depth_dict, probe_profiles_dict


# ---------------------------------------------------------------------------
# Parallel execution wrapper
# ---------------------------------------------------------------------------

def _run_parallel(
    bam_files: list[str],
    args,
    regions_dict: dict,
    plot_intervals_cache: dict | None,
    output_dir: str,
    coverage_dfs: list,
    breakpoints_data: list,
    pooled_depth: dict,
    probe_profiles_all: dict,
) -> None:
    """Distribute per-BAM work across a ``ProcessPoolExecutor``."""
    work_items = [
        (bam_file, args, regions_dict, plot_intervals_cache, output_dir)
        for bam_file in bam_files
    ]

    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(_process_one_bam, item): item[0]
            for item in work_items
        }
        for future in futures:  # as_completed preserves error propagation
            bam_file = futures[future]
            try:
                bam_name, cov_df, bkps, pd_dict, pp_dict = future.result()
                if cov_df is not None:
                    coverage_dfs.append(cov_df)
                breakpoints_data.extend(bkps)
                for locus, entries in pd_dict.items():
                    pooled_depth.setdefault(locus, []).extend(entries)
                for locus, sample_profiles in pp_dict.items():
                    probe_profiles_all.setdefault(locus, {}).update(sample_profiles)
            except Exception as exc:
                print(f"ERROR processing {bam_file}: {exc}")


# ---------------------------------------------------------------------------
# Output helpers
# ---------------------------------------------------------------------------

def _format_row(row) -> str:
    """Safely format one row of coverage statistics for the TSV output."""
    def _f(val, digits=2):
        try:
            return str(round(float(val), digits))
        except (ValueError, TypeError):
            return "NA"

    mean = _f(row["mean"])
    sd = _f(row["sd"])
    median = _f(row["median"])
    try:
        covered = str(int(row["coveredbases"])) if row["coveredbases"] != "NA" else "NA"
    except (ValueError, TypeError):
        covered = "NA"
    # RCN keeps 4 decimals: one decimal hides differences such as 0.04 against 0.08 (e.g. heteroplasmic deletions)
    norm = _f(row["normalised"], 4)
    cell = f"{mean};{sd};{median};{covered};{norm}"
    if "rcn_low" in row.index:
        cell += f";{_f(row['rcn_low'], 4)};{_f(row['rcn_high'], 4)}"
    return cell


# ---------------------------------------------------------------------------
# Normalisation locus and GC check
# ---------------------------------------------------------------------------

NORM_RANDOM = "norm_random"


def _set_normalisation(args, regions_dict, output_dir: str, out_name: str) -> str:
    """Add the normalisation regions to *regions_dict*, point ``args.norm`` at them, check GC.

    Returns a one-line description of the normalisation for the coverage header.
    """
    if args.norm_random is not None:
        if not args.fasta:
            sys.exit("--norm-random needs --fasta to draw the windows.")
        if NORM_RANDOM in regions_dict:
            sys.exit(f"The region file already has a locus named {NORM_RANDOM}.")
        unknown = [name for name in args.norm_ignore if name not in regions_dict]
        if unknown:
            sys.exit(f"--norm-ignore loci not in the region file: {', '.join(unknown)}")
        target_gc = None
        if args.norm_gc == "match":
            targets = {k: v for k, v in regions_dict.items() if k not in args.norm_ignore}
            gcs = [g for g, _ in loci_gc(args.fasta, targets).values() if not np.isnan(g)]
            if not gcs:
                sys.exit("No locus left to match in the fasta, so no GC to match.")
            target_gc = float(np.median(gcs))
        exclude = [
            (chrom, s - args.norm_margin, e + args.norm_margin)
            for region_list in regions_dict.values()
            for chrom, s, e in (parse_region(r) for r in region_list if r)
        ]
        if args.norm_exclude:
            exclude += read_bed_intervals(args.norm_exclude)
        windows = draw_norm_windows(
            args.fasta, args.norm_random, args.norm_size, exclude,
            target_gc=target_gc, chroms=args.norm_chroms, seed=args.norm_seed,
        )
        if not windows:
            sys.exit("No normalisation window could be drawn.")
        if len(windows) < args.norm_random:
            print(f"Warning: only {len(windows)} of {args.norm_random} normalisation windows "
                  "matched; the genome has few windows of this GC.")
        bed = os.path.join(output_dir, f"{out_name}_norm_windows.bed")
        write_windows_bed(bed, windows, NORM_RANDOM)
        print(f"Normalisation windows saved to {bed}")
        regions_dict[NORM_RANDOM] = [f"{c}:{s + 1}-{e}" for c, s, e in windows]
        args.norm = NORM_RANDOM
        ignored = f" but {', '.join(args.norm_ignore)}" if args.norm_ignore else ""
        gc_note = (f", GC {target_gc:.3f} +/- 0.02 (median of the loci{ignored})"
                   if target_gc is not None else "")
        note = (f"{len(windows)} random windows of {args.norm_size} bp{gc_note}, "
                f"seed {args.norm_seed}, in {os.path.basename(bed)}")
    elif os.path.isfile(args.norm):
        norm_regions = get_regions(args.norm)
        if not norm_regions:
            sys.exit(f"No region in the normalisation file {args.norm}.")
        name = next(iter(norm_regions))
        if name in regions_dict:
            sys.exit(f"The normalisation file and the region file share the locus name {name}.")
        regions_dict[name] = [r for region_list in norm_regions.values() for r in region_list]
        note = f"{len(regions_dict[name])} regions from {args.norm}, pooled as {name}"
        args.norm = name
    elif args.norm in regions_dict:
        note = f"locus {args.norm} of {args.region}"
    else:
        raise ValueError(
            f"Normalisation region '{args.norm}' not found in the regions file."
        )

    if args.fasta:
        gc_by_locus = loci_gc(args.fasta, regions_dict)
        if np.isnan(gc_by_locus[args.norm][0]):
            print("Warning: the normalisation regions are not in the fasta; GC check skipped.")
            return note
        gc_out = os.path.join(output_dir, f"{out_name}_gc.tsv")
        flagged = write_gc_report(gc_out, gc_by_locus, args.norm, args.gc_warn)
        print(f"GC report saved to {gc_out}")
        for name in flagged:
            print(f"Warning: GC of {name} differs from {args.norm} by more than "
                  f"{args.gc_warn:g} points (regions +/- {GC_PAD} bp); its RCN may carry a GC bias.")
    else:
        print("No --fasta given: GC of the loci not checked.")
    return note


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run(args) -> None:
    """Entry point for ``ardu coverage``."""
    start_time = time.time()
    configure_read_filters(args.min_mapq, args.min_baseq)
    date_str = time.strftime("%Y-%m-%d")
    ArDuBase = os.path.splitext(os.path.basename(args.outfile))[0]
    ArDuBaseRun = f"{ArDuBase}-{date_str}"
    # the run folder goes next to the -o prefix, so -o can point to another directory
    out_parent = os.path.dirname(args.outfile)
    out_name = os.path.basename(args.outfile)
    if out_parent and not os.path.isdir(out_parent):
        sys.exit(f"Output directory does not exist: {out_parent}")
    # a rerun with the same prefix on the same day gets a numbered suffix instead of overwriting
    output_dir = os.path.join(out_parent, f"ArDuRun-{ArDuBaseRun}")
    k = 2
    while os.path.exists(output_dir):
        output_dir = os.path.join(out_parent, f"ArDuRun-{ArDuBaseRun}_{k}")
        k += 1
    output_dir += "/"
    os.makedirs(output_dir)
    print(f"Output directory: {output_dir}")

    with open(args.bam) as fh:
        bam_files = fh.read().splitlines()

    regions_dict = get_regions(args.region)
    norm_note = _set_normalisation(args, regions_dict, output_dir, out_name)

    plot_intervals_cache = None
    if getattr(args, "plot_interval", None):
        plot_intervals_cache = parse_plot_intervals(args.plot_interval)

    coverage_dfs: list[pd.DataFrame] = []
    breakpoints_data: list[dict] = []
    pooled_depth: dict = {}
    probe_profiles_all: dict = {}

    bam_names_without_ext = [
        os.path.basename(b).removesuffix(".bam") for b in bam_files
    ]

    # ------------------------------------------------------------------
    # Per-BAM processing (serial or parallel)
    # ------------------------------------------------------------------
    if getattr(args, "threads", 1) > 1:
        _run_parallel(
            bam_files, args, regions_dict, plot_intervals_cache, output_dir,
            coverage_dfs, breakpoints_data, pooled_depth, probe_profiles_all,
        )
    else:
        for bam_file in bam_files:
            try:
                _, cov_df, bkps, pd_dict, pp_dict = _process_one_bam(
                    (bam_file, args, regions_dict, plot_intervals_cache, output_dir)
                )
                if cov_df is not None:
                    coverage_dfs.append(cov_df)
                breakpoints_data.extend(bkps)
                for locus, entries in pd_dict.items():
                    pooled_depth.setdefault(locus, []).extend(entries)
                for locus, sample_profiles in pp_dict.items():
                    probe_profiles_all.setdefault(locus, {}).update(sample_profiles)
            except Exception as exc:
                print(f"{bam_file}: {exc}")

    # ------------------------------------------------------------------
    # Pooled plots (after all BAMs processed)
    # ------------------------------------------------------------------
    if args.plot and getattr(args, "plot_pooled", False) and pooled_depth:
        print("Generating pooled plots...")
        for locus, sample_list in pooled_depth.items():
            make_pooled_plot(
                locus, sample_list, args, output_dir, breakpoints_data,
                probe_profiles=(
                    probe_profiles_all.get(locus)
                    if getattr(args, "probe_profile", False)
                    else None
                ),
            )
        if getattr(args, "probe_profile", False) and getattr(args, "probe_profile_out", False):
            for locus, profiles in probe_profiles_all.items():
                totalspan = get_totalspan(regions_dict, locus)
                if totalspan is None:
                    continue
                chrom_p, start_p, end_p = totalspan
                first_bam = bam_files[0]
                pos_arr, _ = build_probe_profile(
                    first_bam, chrom_p, start_p, end_p,
                    ref_median=1.0,  # positions only; ratios ignored here
                    probe_size=args.probe_size,
                    probe_spacing=args.probe_spacing,
                    n_probes_each_side=args.probe_n_each_side,
                )
                from ardu.probe import cluster_probe_profiles
                labels = cluster_probe_profiles(
                    profiles,
                    threshold=args.probe_threshold,
                    distance_threshold=args.probe_cluster_dist,
                )
                write_probe_profiles_tsv(
                    profiles, pos_arr, labels,
                    threshold=args.probe_threshold,
                    output_path=os.path.join(output_dir, f"{locus}_probe_profiles.tsv"),
                )

    # ------------------------------------------------------------------
    # Coverage output
    # ------------------------------------------------------------------
    if coverage_dfs:
        combined_df = pd.concat(coverage_dfs, ignore_index=True)
        combined_df["values"] = combined_df.apply(_format_row, axis=1)
        pivot_df = combined_df.pivot(
            index="locus", columns="bam_file", values="values"
        ).reset_index()
        cols = ["locus"] + [b for b in bam_names_without_ext if b in pivot_df.columns]
        pivot_df = pivot_df[cols]

        cov_out = os.path.join(output_dir, f"{out_name}_coverage.tsv")
        with open(cov_out, "w") as fh:
            fh.write("# ArDu coverage outfile\n")
            fh.write(f"# Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
            if args.rcn_ci > 0:
                fh.write("# Format: uDoC;sdDoC;medDoC;CovBases;RCN;RCN_low;RCN_high\n")
            else:
                fh.write("# Format: uDoC;sdDoC;medDoC;CovBases;RCN\n")
            fh.write("# uDoC:     mean Depth of Coverage over the entire target span.\n")
            fh.write("# sdDoC:    standard deviation of Depth of Coverage.\n")
            fh.write("# medDoC:   median Depth of Coverage over the entire locus span.\n")
            fh.write("# CovBases: total number of covered bases on locus span.\n")
            fh.write("# RCN:      relative copy number (median depth of the locus / median depth "
                     "of the normalisation locus, i.e. copies per normalisation-locus copy).\n")
            fh.write(f"# Normalisation: {norm_note}\n")
            if args.rcn_ci > 0:
                fh.write(f"# RCN_low, RCN_high: 95% interval of RCN from {args.rcn_ci} block bootstrap "
                         f"replicates ({args.ci_block} bp blocks, seed {args.ci_seed}).\n")
            pivot_df.to_csv(fh, sep="\t", index=False)
        print(f"Coverage data saved to {cov_out}")
    else:
        print("No coverage data to write.")

    # ------------------------------------------------------------------
    # Breakpoints output
    # ------------------------------------------------------------------
    if breakpoints_data:
        bp_df = pd.DataFrame(breakpoints_data)
        bp_out = os.path.join(output_dir, f"{out_name}_breakpoints.tsv")
        bp_df.to_csv(bp_out, sep="\t", index=False)
        print(f"Breakpoints data saved to {bp_out}")

    # ------------------------------------------------------------------
    # Mutation genotyping output
    # ------------------------------------------------------------------
    if getattr(args, "mutation", None):
        mutation_data: dict = {}
        with open(args.mutation) as mf:
            for line in mf:
                parts = line.strip().split("\t")
                chromosome = parts[0]
                position = parts[1]
                mutation_name = parts[2] if len(parts) > 2 else "noname"
                mutation_id = f"{chromosome}:{position}:{mutation_name}"
                mutation_data[mutation_id] = {}
                for bam_file in bam_files:
                    try:
                        nucleotide_counts, total_depth = calculate_nucleotide_counts(
                            bam_file, chromosome, position
                        )
                        mutation_data[mutation_id][bam_file] = {
                            "A": nucleotide_counts["A"],
                            "T": nucleotide_counts["T"],
                            "C": nucleotide_counts["C"],
                            "G": nucleotide_counts["G"],
                            "depth": total_depth,
                        }
                        if args.mutation_indels:
                            mutation_data[mutation_id][bam_file].update(
                                calculate_indel_counts(bam_file, chromosome, position)
                            )
                    except Exception as exc:
                        print(f"Error processing {mutation_id} in {bam_file}: {exc}")
                        mutation_data[mutation_id][bam_file] = {
                            "A": "NA", "T": "NA", "C": "NA", "G": "NA", "depth": "NA"
                        }

        mut_out = os.path.join(output_dir, f"{out_name}_mutations.tsv")
        with open(mut_out, "w") as out_file:
            out_file.write("mutation_id\t" + "\t".join(bam_files) + "\n")
            for mutation_id, bam_counts in mutation_data.items():
                row = [mutation_id]
                for bam_file in bam_files:
                    c = bam_counts.get(bam_file, {})
                    cell = (
                        f"A={c.get('A','NA')};T={c.get('T','NA')};"
                        f"C={c.get('C','NA')};G={c.get('G','NA')};"
                        f"depth={c.get('depth','NA')}"
                    )
                    if args.mutation_indels:
                        runs = c.get("runs") or {}
                        run_str = ",".join(
                            f"{k}:{v}" for k, v in sorted(runs.items(), key=lambda kv: (kv[0][0], int(kv[0][1:])))
                        ) or "NA"
                        cell += f";ins={c.get('ins','NA')};del={c.get('del','NA')};runs={run_str}"
                    row.append(cell)
                out_file.write("\t".join(row) + "\n")
        print(f"Mutation data saved to {mut_out}")

    elapsed = round(time.time() - start_time, 1)
    print(f"Elapsed time: {elapsed} seconds")
