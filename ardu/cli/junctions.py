"""
``ardu junctions`` — extract softclips, find junction reads via Aho-Corasick,
detect bridging pairs, run BLAST, and produce insert-size / softclip plots.

Corresponds to JunctionSequenceBySoftclip-Finder.py, refactored as a
subcommand with a ``register()`` / ``run()`` interface.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
import pysam

matplotlib.use("Agg")

from ardu.junctionlib import (
    build_softclip_automaton,
    find_reads_with_multiple_softclips_aho,
    write_junction_reads_fasta,
    write_junction_reads_fastq,
)


# ---------------------------------------------------------------------------
# Argparse registration
# ---------------------------------------------------------------------------

def register(subparsers) -> None:
    """Add the ``junctions`` sub-command to *subparsers*."""
    parser = subparsers.add_parser(
        "junctions",
        help=(
            "Softclip extraction, junction-read detection, bridging pairs, "
            "BLAST, insert-size and softclip plots."
        ),
        description=(
            "ArDu junctions — full junction-evidence pipeline: softclip "
            "extraction, Aho-Corasick junction search, bridging read detection, "
            "optional BLAST, and diagnostic plots."
        ),
    )

    parser.add_argument("-b", "--bam_list", required=True,
                        help="File listing BAM paths, one per line.")
    parser.add_argument("-i", "--input", required=True,
                        help="Breakpoints TSV as produced by ``ardu coverage``.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output file prefix.")
    parser.add_argument("-s", "--size", type=int, default=30,
                        help="Minimum soft-clip length to extract. Default=30.")
    parser.add_argument("-e", "--extension", type=int, default=30,
                        help="bp to extend around each breakpoint position. Default=30.")
    parser.add_argument("--blast", metavar="REFERENCE",
                        help="Reference genome FASTA; triggers BLAST run + filtering.")
    parser.add_argument("--pairs", nargs="+", metavar="N-M",
                        help="Allowed breakpoint pairs, e.g. --pairs 1-2 3-4")
    parser.add_argument("--junction", action="store_true",
                        help="Run Aho-Corasick junction-read search.")
    parser.add_argument("--junction-fastq", action="store_true",
                        help="Write junction reads as FASTQ (default: FASTA).")
    parser.add_argument("--plot-insert", action="store_true",
                        help="Plot insert-size distributions.")
    parser.add_argument("--plot-softclip", action="store_true",
                        help="Plot softclip position distributions.")
    parser.add_argument("--plot-mode",
                        choices=["per-bam", "pooled", "both"],
                        default="per-bam",
                        help="Plotting mode. Default=per-bam.")

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def _extract_breakpoint_from_qseqid(qseqid: str) -> int | None:
    match = re.search(r"_(\d+)_", qseqid)
    return int(match.group(1)) if match else None


def _parse_allowed_pairs(pairs_arg: list[str] | None) -> set[tuple[int, int]] | None:
    """Parse ``["1-2", "3-4"]`` into a set of both-direction int tuples."""
    if not pairs_arg:
        return None
    allowed: set[tuple[int, int]] = set()
    for token in pairs_arg:
        token = token.strip()
        if not token or "-" not in token:
            raise ValueError(f"Invalid --pairs value '{token}'. Use e.g. --pairs 1-2 3-4")
        i, j = map(int, token.split("-", 1))
        allowed.add((i, j))
        allowed.add((j, i))
    return allowed


def _canonical_pair_key(a: int, b: int) -> tuple[int, int]:
    return tuple(sorted((a, b)))  # type: ignore[return-value]


def _filter_blast_hits(
    blast_output: str,
    breakpoints_df: pd.DataFrame,
    extension: int,
    filtered_output: str,
    allowed_pairs: set[tuple[int, int]] | None,
) -> None:
    regions_by_bkp: dict[int, tuple] = {}
    for _, row in breakpoints_df.iterrows():
        regions_by_bkp[row["breakpoint_number"]] = (
            row["chromosome"],
            row["position"] - extension,
            row["position"] + extension,
        )

    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast_df = pd.read_csv(blast_output, sep="\t", names=cols)

    filtered_hits = []
    for _, row in blast_df.iterrows():
        source_bkp = _extract_breakpoint_from_qseqid(str(row["qseqid"]))
        if source_bkp is None:
            continue
        sstart = min(row["sstart"], row["send"])
        send = max(row["sstart"], row["send"])
        sseqid = str(row["sseqid"])
        for target_bkp, (chrom, start, end) in regions_by_bkp.items():
            if target_bkp == source_bkp:
                continue
            if allowed_pairs is not None and (source_bkp, target_bkp) not in allowed_pairs:
                continue
            if sseqid == chrom and not (send < start or sstart > end):
                filtered_hits.append(row)
                break

    pd.DataFrame(filtered_hits).to_csv(filtered_output, sep="\t", index=False, header=False)
    print(f"Filtered BLAST results saved to {filtered_output}")


def _run_blast(query_fasta: str, reference: str, blast_output: str) -> None:
    if not os.path.exists(f"{reference}.nhr"):
        print("Creating BLAST database...")
        subprocess.run(["makeblastdb", "-in", reference, "-dbtype", "nucl"], check=True)
    print("Running BLAST...")
    subprocess.run(
        ["blastn", "-query", query_fasta, "-db", reference,
         "-out", blast_output, "-outfmt", "6"],
        check=True,
    )
    print(f"BLAST results saved in {blast_output}")


def _extract_softclips(
    bam_file: str,
    region: tuple,
    min_clip: int,
    fasta_output: str,
    tsv_output: str,
    breakpoint: int,
    inserts_by_bkp: dict,
    regions_by_bkp: dict,
    bridges_writer,
    softclip_positions_by_bkp: dict,
) -> dict[int, list[str]]:
    """Extract softclips and return ``{bkp: [clip_sequence, ...]}``.

    Also writes to the shared FASTA / TSV outputs, appends insert sizes,
    records softclip positions, and logs bridging read pairs.
    """
    softclips_by_bkp: dict[int, list[str]] = {breakpoint: []}

    def _has_valid_cigar(read) -> bool:
        return read.cigar is not None and len(read.cigar) > 0

    with pysam.AlignmentFile(bam_file, "rb") as bam, \
         open(fasta_output, "a") as fasta_out, \
         open(tsv_output, "a") as tsv_out:

        for read in bam.fetch(region[0], region[1], region[2]):
            if read.query_sequence is None or not _has_valid_cigar(read):
                continue

            mapping_position = read.reference_start + 1
            chrom = read.reference_name
            cigar = read.cigarstring
            mate_position = (
                read.next_reference_start + 1
                if not read.mate_is_unmapped
                else None
            )
            insert_size = (
                abs(mapping_position - mate_position)
                if mate_position is not None
                and read.reference_name == read.next_reference_name
                else None
            )

            inserts_by_bkp.setdefault(breakpoint, [])
            if insert_size is not None:
                inserts_by_bkp[breakpoint].append(insert_size)

            # LEFT softclip
            if read.cigar[0][0] == 4 and read.cigar[0][1] >= min_clip:
                left_clip = read.query_sequence[: read.cigar[0][1]]
                softclips_by_bkp[breakpoint].append(left_clip)
                fasta_out.write(
                    f">{chrom}:{mapping_position}_{read.query_name}_{breakpoint}_L\n"
                    f"{left_clip}\n"
                )
                tsv_out.write(
                    f"{bam_file}\t{read.query_name}\t{mapping_position}\t"
                    f"{mapping_position}\tleftmost_clip\t{cigar}\t{mate_position}\t"
                    f"{insert_size}\t{region}\t{breakpoint}\n"
                )
                softclip_positions_by_bkp.setdefault(breakpoint, []).append(
                    mapping_position
                )

            # RIGHT softclip
            if read.cigar[-1][0] == 4 and read.cigar[-1][1] >= min_clip:
                right_clip = read.query_sequence[-read.cigar[-1][1]:]
                clip_position = read.reference_end
                softclips_by_bkp[breakpoint].append(right_clip)
                fasta_out.write(
                    f">{chrom}:{clip_position}_{read.query_name}_{breakpoint}_R\n"
                    f"{right_clip}\n"
                )
                tsv_out.write(
                    f"{bam_file}\t{read.query_name}\t{mapping_position}\t"
                    f"{clip_position}\trightmost_clip\t{cigar}\t{mate_position}\t"
                    f"{insert_size}\t{region}\t{breakpoint}\n"
                )
                softclip_positions_by_bkp.setdefault(breakpoint, []).append(
                    clip_position
                )

            # Bridging read detection
            if (
                bridges_writer
                and mate_position is not None
                and read.reference_name == read.next_reference_name
            ):
                for target_bkp, (t_chrom, t_start, t_end) in regions_by_bkp.items():
                    if target_bkp == breakpoint:
                        continue
                    if t_start <= mate_position <= t_end:
                        bridges_writer.write(
                            f"{read.query_name}\t{read.query_name}\t"
                            f"{chrom}:{mapping_position}\t{chrom}:{mate_position}\t"
                            f"{insert_size}\t{breakpoint}-{target_bkp}\n"
                        )
                        break

    return softclips_by_bkp


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _plot_insert_sizes(
    inserts_by_pair: dict,
    output_path: str,
    bkp_to_expected_pos: dict | None = None,
    library_insert: int = 500,
) -> None:
    """Plot insert-size distributions per breakpoint pair.

    When *bkp_to_expected_pos* is provided, the expected bridging insert size
    is computed from the distance between the two breakpoints.  Inserts beyond
    (duplication_size + library_insert + 10% of duplication_size) are treated
    as likely misalignments: they are excluded from the plot but counted and
    reported in the panel title.  A vertical line marks the expected insert.
    """
    num = max(1, len(inserts_by_pair))
    cols = min(num, 3)
    rows = math.ceil(num / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows))
    axes = np.array(axes).reshape(-1)

    for ax, pair_key in zip(axes, inserts_by_pair):
        data = inserts_by_pair[pair_key]
        if not data:
            ax.axis("off")
            continue

        # Compute expected insert and plausibility cap from breakpoint positions
        expected_insert = None
        max_plausible = None
        if bkp_to_expected_pos and len(pair_key) == 2:
            b1, b2 = pair_key
            pos1 = bkp_to_expected_pos.get(b1)
            pos2 = bkp_to_expected_pos.get(b2)
            if pos1 and pos2:
                dup_size = abs(pos2 - pos1)
                tolerance = dup_size * 0.10
                expected_insert = dup_size + library_insert
                max_plausible = expected_insert + tolerance

        # Split data into plausible and outliers
        if max_plausible is not None:
            plausible = [x for x in data if x <= max_plausible]
            outliers = [x for x in data if x > max_plausible]
        else:
            plausible = data
            outliers = []

        if not plausible:
            ax.axis("off")
            continue

        ax.hist(
            plausible,
            bins=min(100, max(10, int(np.sqrt(len(plausible))))),
            edgecolor="black",
            color="steelblue",
            label=f"n={len(plausible)}",
        )
        ax.set_yscale("log")

        if expected_insert is not None:
            ax.axvline(
                expected_insert, color="red", linestyle="--", linewidth=1.5,
                label=f"expected: {expected_insert:,.0f} bp",
            )

        pair_label = "-".join(map(str, pair_key))
        title = f"Pair {pair_label}"
        if outliers:
            title += f" ({len(outliers)} outliers excluded)"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Insert size (bp)", fontsize=9)
        ax.set_ylabel("Count (log)", fontsize=9)
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}")
        )
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
        ax.legend(fontsize=8)

    for ax in axes[len(inserts_by_pair):]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path}")


def _plot_softclip_positions(
    softclip_positions_by_bkp: dict,
    bkp_to_expected_pos: dict,
    output_path: str,
    bam_name: str = "",
    pair_keys: list | None = None,
) -> None:
    """Plot softclip pileup per breakpoint as bar charts with actual genomic coordinates.

    One panel per breakpoint (per-BAM mode) or per pair (pooled mode).
    Each panel shows:
      - bar chart of softclip counts at each position
      - blue dashed line: mode (position with most softclips)
      - red line: predicted breakpoint position from coverage analysis
      - offset between mode and predicted shown in the title
    """
    # Build per-panel data: list of (panel_title, {bkp: positions_list})
    if pair_keys is not None:
        # pooled mode: one panel per pair, positions from both breakpoints overlaid
        data_items = []
        for pair_key in pair_keys:
            panel_data = {
                bkp: softclip_positions_by_bkp.get(bkp, [])
                for bkp in pair_key
            }
            data_items.append((pair_key, panel_data))
    else:
        # per-BAM mode: one panel per breakpoint
        data_items = [
            ((bkp,), {bkp: positions})
            for bkp, positions in softclip_positions_by_bkp.items()
        ]

    num = max(1, len(data_items))
    cols = min(num, 3)
    rows = math.ceil(num / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows))
    axes = np.array(axes).reshape(-1)

    colors = ["steelblue", "darkorange", "seagreen", "mediumpurple"]

    for ax, (pair_key, panel_data) in zip(axes, data_items):
        has_data = any(len(v) > 0 for v in panel_data.values())
        if not has_data:
            ax.axis("off")
            continue

        for color_idx, (bkp, positions) in enumerate(panel_data.items()):
            if not positions:
                continue
            counts = Counter(positions)
            xs = sorted(counts)
            ys = [counts[x] for x in xs]
            color = colors[color_idx % len(colors)]

            # Bar chart — width scaled to the density of positions
            span = max(xs) - min(xs) if len(xs) > 1 else 1
            bar_width = max(1, span // max(len(xs), 1))
            ax.bar(xs, ys, width=bar_width, color=color, alpha=0.7,
                   label=f"BKP {bkp} (n={len(positions)})")

            # Mode line
            mode = max(counts, key=counts.__getitem__)
            ax.axvline(mode, color=color, linestyle="--", linewidth=1.5,
                       label=f"BKP {bkp} mode: {mode:,}")

            # Predicted breakpoint line + offset annotation
            exp = bkp_to_expected_pos.get(bkp)
            if exp:
                ax.axvline(exp, color="red", linestyle="-", linewidth=1.2,
                           label=f"BKP {bkp} predicted: {exp:,}")
                offset = mode - exp
                ax.annotate(
                    f"offset: {offset:+,} bp",
                    xy=(exp, ax.get_ylim()[1]),
                    xytext=(5, -12),
                    textcoords="offset points",
                    color="red", fontsize=8,
                )

        chrom = ""
        for bkp in pair_key:
            exp = bkp_to_expected_pos.get(bkp)
            if exp:
                # get chrom from regions if available (best effort)
                break
        bkp_label = "-".join(map(str, pair_key))
        title = f"BKP {bkp_label}"
        if bam_name:
            title = f"{bam_name} | {title}"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Genomic position", fontsize=9)
        ax.set_ylabel("Softclip count", fontsize=9)
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}")
        )
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
        ax.legend(fontsize=7, loc="upper right")

    for ax in axes[num:]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path}")


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run(args) -> None:
    """Entry point for ``ardu junctions``."""
    fasta_output = f"{args.output}_soft_clipped_sequences.fasta"
    tsv_output = f"{args.output}.tsv"
    blast_output = f"{args.output}_BLAST_results.txt"
    bridges_output = f"{args.output}_bridging_pairs.tsv"

    for path in (fasta_output, tsv_output):
        if os.path.exists(path):
            os.remove(path)

    allowed_pairs = _parse_allowed_pairs(getattr(args, "pairs", None))
    canonical_pairs: list[tuple[int, int]] = []
    if allowed_pairs:
        seen: set = set()
        for a, b in allowed_pairs:
            key = _canonical_pair_key(a, b)
            if key not in seen:
                canonical_pairs.append(key)
                seen.add(key)

    bam_mapping: dict[str, str] = {}
    with open(args.bam_list) as fh:
        for line in fh:
            full = line.strip()
            if not full:
                continue
            base = os.path.basename(full)
            if base.endswith(".bam"):
                bam_mapping[base] = full

    # Accept either a direct path to a breakpoints TSV or an ArDuRun directory
    import glob as _glob
    if os.path.isdir(args.input):
        candidates = _glob.glob(os.path.join(args.input, "*_breakpoints.tsv"))
        if len(candidates) == 0:
            raise FileNotFoundError(
                f"No *_breakpoints.tsv found in {args.input}. "
                "Run 'ardu coverage' with --breakpoint first."
            )
        if len(candidates) > 1:
            raise ValueError(
                f"Multiple breakpoints files found in {args.input}:\n"
                + "\n".join(candidates)
                + "\nPlease specify one directly with -i."
            )
        breakpoints_path = candidates[0]
        print(f"Using breakpoints file: {breakpoints_path}")
    else:
        breakpoints_path = args.input

    breakpoints_df = pd.read_csv(breakpoints_path, sep="\t")
    breakpoints_df["bam_file"] = breakpoints_df["bam_file"].apply(
        lambda x: os.path.basename(x) + ".bam"
    )

    regions_by_bkp = {
        row["breakpoint_number"]: (
            row["chromosome"],
            row["position"] - args.extension,
            row["position"] + args.extension,
        )
        for _, row in breakpoints_df.iterrows()
    }
    bkp_to_expected_pos = {
        row["breakpoint_number"]: row["position"]
        for _, row in breakpoints_df.iterrows()
    }

    with open(tsv_output, "w") as out:
        out.write(
            "bam_file\tRead_name\tPosition\tClipping_position\t"
            "Clipping_end\tCIGAR\tMate_position\tInsert_size\tRegion\tBreakpoint\n"
        )
    with open(bridges_output, "w") as bf:
        bf.write("read_name\tmate_name\tread_pos\tmate_pos\tinsert_size\tbreakpoint_bridge\n")

    pooled_inserts_by_pair: dict = {}
    pooled_softclip_by_pair: dict = {}

    # ------------------------------------------------------------------
    # Per-BAM processing
    # ------------------------------------------------------------------
    for bam_basename, bam_file in bam_mapping.items():
        print(f"\n=== PROCESSING BAM: {bam_basename} ===")
        bam_name = bam_basename.removesuffix(".bam")

        softclips_by_bkp: dict[int, list[str]] = {}
        softclip_positions_by_bkp: dict[int, list[int]] = {}
        inserts_by_bkp: dict[int, list[int]] = {}

        bkp_subset = breakpoints_df[breakpoints_df["bam_file"] == bam_basename]

        for _, row in bkp_subset.iterrows():
            bkp = row["breakpoint_number"]
            pos = row["position"]
            region = (row["chromosome"], pos - args.extension, pos + args.extension)

            with open(bridges_output, "a") as bridges_fh:
                bkp_clips = _extract_softclips(
                    bam_file, region, args.size,
                    fasta_output, tsv_output, bkp,
                    inserts_by_bkp, regions_by_bkp,
                    bridges_fh, softclip_positions_by_bkp,
                )
            softclips_by_bkp.setdefault(bkp, []).extend(bkp_clips.get(bkp, []))

        # Pooled accumulation
        if args.plot_mode in ("pooled", "both"):
            if canonical_pairs:
                for key in canonical_pairs:
                    b1, b2 = key
                    pooled_inserts_by_pair.setdefault(key, [])
                    pooled_inserts_by_pair[key].extend(inserts_by_bkp.get(b1, []))
                    pooled_inserts_by_pair[key].extend(inserts_by_bkp.get(b2, []))
                    pooled_softclip_by_pair.setdefault(key, [])
                    pooled_softclip_by_pair[key].extend(softclip_positions_by_bkp.get(b1, []))
                    pooled_softclip_by_pair[key].extend(softclip_positions_by_bkp.get(b2, []))
            else:
                for bkp, sizes in inserts_by_bkp.items():
                    pooled_inserts_by_pair.setdefault((bkp,), []).extend(sizes)
                for bkp, poss in softclip_positions_by_bkp.items():
                    pooled_softclip_by_pair.setdefault((bkp,), []).extend(poss)

        # ---- Junction search ----
        # The search region is derived directly from the breakpoint windows
        # already stored in regions_by_bkp — no need for a separate argument.
        if args.junction:
            print(f"Building automaton for {bam_name}...")
            allowed_bkps = (
                {b for pair in allowed_pairs for b in pair} if allowed_pairs else None
            )
            automaton = build_softclip_automaton(softclips_by_bkp, allowed_bkps)

            # Use the union span of all relevant breakpoint regions as the
            # fetch window so we catch reads spanning between breakpoints.
            relevant_bkps = list(
                (allowed_bkps or set(regions_by_bkp.keys()))
                & set(regions_by_bkp.keys())
            )
            if relevant_bkps:
                chroms = [regions_by_bkp[b][0] for b in relevant_bkps]
                starts = [regions_by_bkp[b][1] for b in relevant_bkps]
                ends   = [regions_by_bkp[b][2] for b in relevant_bkps]
                # Only fetch as a single span when all breakpoints are on the
                # same chromosome; otherwise scan each region individually.
                if len(set(chroms)) == 1:
                    search_regions = [(chroms[0], min(starts), max(ends))]
                else:
                    search_regions = list(regions_by_bkp.values())
            else:
                search_regions = list(regions_by_bkp.values())

            all_matches = []
            for search_region in search_regions:
                all_matches.extend(
                    find_reads_with_multiple_softclips_aho(
                        bam_file, automaton, search_region, allowed_pairs
                    )
                )

            # Deduplicate by read name + breakpoint pair
            seen_keys: set = set()
            matches = []
            for m in all_matches:
                key = (m["read_name"], m["breakpoints"])
                if key not in seen_keys:
                    seen_keys.add(key)
                    matches.append(m)

            output_dir = os.path.dirname(args.output) or "."
            if getattr(args, "junction_fastq", False):
                write_junction_reads_fastq(bam_file, matches, output_dir, bam_name)
            else:
                write_junction_reads_fasta(matches, output_dir, bam_name)

        # ---- Per-BAM plots ----
        if args.plot_insert and args.plot_mode in ("per-bam", "both"):
            if canonical_pairs:
                inserts_by_pair = {
                    key: inserts_by_bkp.get(key[0], []) + inserts_by_bkp.get(key[1], [])
                    for key in canonical_pairs
                }
            else:
                inserts_by_pair = {(bkp,): ins for bkp, ins in inserts_by_bkp.items()}
            _plot_insert_sizes(
                inserts_by_pair,
                f"{args.output}_{bam_name}_insert_dist.png",
                bkp_to_expected_pos=bkp_to_expected_pos,
            )

        if args.plot_softclip and args.plot_mode in ("per-bam", "both"):
            _plot_softclip_positions(
                softclip_positions_by_bkp, bkp_to_expected_pos,
                f"{args.output}_{bam_name}_softclip_dist.png",
                bam_name=bam_name,
            )

    # ------------------------------------------------------------------
    # BLAST
    # ------------------------------------------------------------------
    if args.blast:
        print("\n=== RUNNING BLAST ON ALL SOFTCLIPS ===")
        _run_blast(fasta_output, args.blast, blast_output)
        filtered_output = f"{args.output}_BLAST_filtered_pairs.txt"
        _filter_blast_hits(blast_output, breakpoints_df, args.extension,
                           filtered_output, allowed_pairs)

    # ------------------------------------------------------------------
    # Pooled plots
    # ------------------------------------------------------------------
    if args.plot_mode in ("pooled", "both"):
        if args.plot_insert and pooled_inserts_by_pair:
            print("Plotting pooled insert-size distributions by pair...")
            _plot_insert_sizes(
                pooled_inserts_by_pair,
                f"{args.output}_pooled_insert_dist_by_pair.png",
                bkp_to_expected_pos=bkp_to_expected_pos,
            )
        if args.plot_softclip and pooled_softclip_by_pair:
            print("Plotting pooled softclip distributions by pair...")
            _plot_softclip_positions(
                pooled_softclip_by_pair, bkp_to_expected_pos,
                f"{args.output}_pooled_softclip_dist_by_pair.png",
                bam_name="pooled",
                pair_keys=list(pooled_softclip_by_pair.keys()),
            )

    print(f"\nBridging pairs saved to {bridges_output}")
    print("Done.")
