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

from ardu.bam import read_pair_info
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
    parser.add_argument("--sv-type", choices=["auto", "deletion", "duplication"],
                        default="auto",
                        help="SV type used for the expected insert of bridging pairs. "
                             "auto picks the type the pair orientations support "
                             "(FR: deletion, RF: tandem duplication), which assumes "
                             "the SV is longer than the library insert. Default=auto.")
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
    pairs: dict,
    softclip_positions_by_bkp: dict,
) -> dict[int, list[str]]:
    """Extract softclips and return ``{bkp: [clip_sequence, ...]}``.

    Also writes to the shared FASTA / TSV outputs, records softclip positions,
    and stores in *pairs* ``{read_name: (read_start, PairInfo)}`` for every
    read pair with a mate in the window.
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
            pair = read_pair_info(read)
            if pair is not None:
                pairs[read.query_name] = (read.reference_start, pair)
            insert_size = pair.insert if pair is not None else None

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

    return softclips_by_bkp


# ---------------------------------------------------------------------------
# Insert sizes
# ---------------------------------------------------------------------------

# Pairs further apart than library median + _DISCORDANT_SD * sd are discordant
_DISCORDANT_SD = 4
_MAD_TO_SD = 1.4826


def _median_mad(values: list[int]) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    a = np.asarray(values, dtype=float)
    median = float(np.median(a))
    return median, float(np.median(np.abs(a - median)))


def _library_inserts(pairs_by_key: dict) -> list[int]:
    """Inserts of the proper FR pairs over all windows, each pair counted once."""
    seen: dict[str, int] = {}
    for pairs in pairs_by_key.values():
        for name, (_, pair) in pairs.items():
            if pair.proper and pair.orientation == "FR":
                seen[name] = pair.insert
    return list(seen.values())


def _is_discordant(pair, library_median: float | None, library_mad: float | None) -> bool:
    if pair.orientation != "FR":
        return True
    if library_median is None:
        return not pair.proper
    return pair.insert > library_median + _DISCORDANT_SD * _MAD_TO_SD * library_mad


def _bridging_pairs(pairs_by_key: dict, regions_by_key: dict,
                    library_median, library_mad) -> dict:
    """Discordant pairs with one mate in a breakpoint window and the other in
    another window of the same locus, as
    ``{(locus, b1, b2): {read_name: (read_start, PairInfo)}}`` with ``b1 < b2``.
    A pair seen from both windows is kept once. Concordant pairs are left out,
    since they bridge any two breakpoints closer than the library insert.
    """
    bridges: dict = {}
    for (locus, bkp), pairs in pairs_by_key.items():
        for (t_locus, t_bkp), (_, t_start, t_end) in regions_by_key.items():
            if t_locus != locus or t_bkp == bkp:
                continue
            key = (locus, *_canonical_pair_key(bkp, t_bkp))
            for name, (read_start, pair) in pairs.items():
                if (pair.mate_end > t_start and pair.mate_start < t_end
                        and _is_discordant(pair, library_median, library_mad)):
                    bridges.setdefault(key, {})[name] = (read_start, pair)
    return bridges


def _finish_panel(panel: dict, library_median, library_mad, sv_type: str) -> None:
    """Set the SV type and expected insert of a bridging *panel*.

    A deletion of size S joins two mates S + L apart (L, the library insert),
    in FR. A tandem duplication joins the end of the copy back to its start, so
    the mates map |L - S| apart in RF. With *sv_type* ``auto``, the type is the
    one the bridging pairs support most, which assumes S > L.
    """
    panel["library_median"] = library_median
    panel["sv_type"] = None
    panel["expected"] = None
    if panel["class"] != "bridging" or library_median is None:
        return
    if sv_type == "auto":
        rf = sum(p.orientation == "RF" for p in panel["pairs"])
        fr = sum(p.orientation == "FR" for p in panel["pairs"])
        sv_type = "duplication" if rf > fr else "deletion" if fr else None
    if sv_type is None:
        return
    panel["sv_type"] = sv_type
    size = panel["size"]
    panel["expected"] = (
        size + library_median if sv_type == "deletion" else abs(library_median - size)
    )


def _summary_rows(bam_name: str, library: list[int], library_mad, panels: list[dict]) -> list[dict]:
    def stats(values):
        a = np.asarray(values, dtype=float)
        return dict(n=len(a), median=np.median(a), q1=np.percentile(a, 25),
                    q3=np.percentile(a, 75))

    base = dict(bam=bam_name, sv_type=None, expected=None, mad=None)
    rows = []
    if library:
        rows.append({**base, "locus": None, "panel": "all", "class": "library",
                     "orientation": "FR", **stats(library), "mad": library_mad})
    for panel in panels:
        info = {**base, "locus": panel["locus"], "panel": panel["panel"],
                "class": panel["class"], "sv_type": panel["sv_type"],
                "expected": panel["expected"]}
        by_orientation: dict[str, list[int]] = {}
        for p in panel["pairs"]:
            by_orientation.setdefault(p.orientation, []).append(p.insert)
        if not by_orientation:
            rows.append({**info, "orientation": None, "n": 0})
        for orientation in sorted(by_orientation):
            rows.append({**info, "orientation": orientation,
                         **stats(by_orientation[orientation])})
    return rows


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

_ORIENTATION_COLOURS = {"FR": "steelblue", "RF": "darkorange", "FF": "seagreen", "RR": "mediumpurple"}


def _plot_insert_sizes(panels: list[dict], output_path: str) -> None:
    """Plot the insert sizes of each panel, stacked by pair orientation.

    Bridging panels (one per breakpoint pair) show the pairs joining the two
    breakpoints, with a red line at the insert expected for the SV type.
    Discordant panels (one per breakpoint) show the discordant pairs with a
    mate in the window. The library insert is given in the title.
    """
    num = max(1, len(panels))
    cols = min(num, 3)
    rows = math.ceil(num / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows))
    axes = np.array(axes).reshape(-1)

    for ax, panel in zip(axes, panels):
        if not panel["pairs"]:
            ax.axis("off")
            continue
        inserts = [p.insert for p in panel["pairs"]]
        bins = np.histogram_bin_edges(
            inserts, bins=min(100, max(10, int(np.sqrt(len(inserts)))))
        )
        orientations = [o for o in _ORIENTATION_COLOURS
                        if any(p.orientation == o for p in panel["pairs"])]
        counts, _, _ = ax.hist(
            [[p.insert for p in panel["pairs"] if p.orientation == o] for o in orientations],
            bins=bins, stacked=True, edgecolor="black",
            color=[_ORIENTATION_COLOURS[o] for o in orientations],
            label=[f"{o} (n={sum(p.orientation == o for p in panel['pairs'])})"
                   for o in orientations],
        )
        # a log axis only helps when counts span orders of magnitude
        log_y = np.max(counts) >= 20
        if log_y:
            ax.set_yscale("log")
        else:
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))

        if panel["expected"] is not None:
            ax.axvline(
                panel["expected"], color="red", linestyle="--", linewidth=1.5,
                label=f"expected {panel['sv_type']}: {panel['expected']:,.0f} bp",
            )

        title = f"{panel['locus']} | BKP {panel['panel']} ({panel['class']} pairs)"
        if panel["library_median"] is not None:
            title += f"\nlibrary insert {panel['library_median']:,.0f} bp"
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Insert size (bp)", fontsize=9)
        ax.set_ylabel("Count (log)" if log_y else "Count", fontsize=9)
        ax.xaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}")
        )
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
        ax.legend(fontsize=8)

    for ax in axes[len(panels):]:
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
    summary_output = f"{args.output}_insert_summary.tsv"

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
        bf.write("bam_file\tlocus\tread_name\tread_pos\tmate_pos\tinsert_size\t"
                 "orientation\tproper_pair\tbreakpoint_bridge\n")

    summary_rows: list[dict] = []
    pooled_library: list[int] = []
    pooled_panels: dict = {}
    pooled_softclip_by_pair: dict = {}

    # ------------------------------------------------------------------
    # Per-BAM processing
    # ------------------------------------------------------------------
    for bam_basename, bam_file in bam_mapping.items():
        print(f"\n=== PROCESSING BAM: {bam_basename} ===")
        bam_name = bam_basename.removesuffix(".bam")

        softclips_by_bkp: dict[int, list[str]] = {}
        softclip_positions_by_bkp: dict[int, list[int]] = {}
        # Breakpoint numbers restart for each sample and locus, so the insert
        # analysis keys its windows by (locus, breakpoint) within this BAM.
        pairs_by_key: dict[tuple, dict] = {}
        regions_by_key: dict[tuple, tuple] = {}
        positions_by_key: dict[tuple, int] = {}

        bkp_subset = breakpoints_df[breakpoints_df["bam_file"] == bam_basename]

        for _, row in bkp_subset.iterrows():
            bkp = row["breakpoint_number"]
            pos = row["position"]
            region = (row["chromosome"], pos - args.extension, pos + args.extension)
            key = (row.get("locus", ""), bkp)
            regions_by_key[key] = region
            positions_by_key[key] = pos

            bkp_clips = _extract_softclips(
                bam_file, region, args.size,
                fasta_output, tsv_output, bkp,
                pairs_by_key.setdefault(key, {}), softclip_positions_by_bkp,
            )
            softclips_by_bkp.setdefault(bkp, []).extend(bkp_clips.get(bkp, []))

        # ---- Insert sizes ----
        library = _library_inserts(pairs_by_key)
        library_median, library_mad = _median_mad(library)
        if library_median is None:
            print(f"No proper FR pairs in the windows of {bam_name}; "
                  "discordant pairs are those not flagged proper by the aligner.")
        elif 2 * args.extension < library_median:
            print(f"Warning: windows of {2 * args.extension} bp are shorter than the "
                  f"library insert ({library_median:,.0f} bp) of {bam_name}; few pairs "
                  "will bridge two breakpoints. Consider a larger -e.")

        bridges = _bridging_pairs(pairs_by_key, regions_by_key, library_median, library_mad)
        with open(bridges_output, "a") as bf:
            for (locus, b1, b2), pairs in bridges.items():
                chrom = regions_by_key[(locus, b1)][0]
                for name, (read_start, pair) in pairs.items():
                    bf.write(
                        f"{bam_basename}\t{locus}\t{name}\t{chrom}:{read_start + 1}\t"
                        f"{chrom}:{pair.mate_start + 1}\t{pair.insert}\t"
                        f"{pair.orientation}\t{pair.proper}\t{b1}-{b2}\n"
                    )

        panels: list[dict] = []
        if canonical_pairs:
            for locus in dict.fromkeys(locus for locus, _ in pairs_by_key):
                for b1, b2 in canonical_pairs:
                    if (locus, b1) not in positions_by_key or (locus, b2) not in positions_by_key:
                        continue
                    panels.append({
                        "locus": locus, "panel": f"{b1}-{b2}", "class": "bridging",
                        "pairs": [p for _, p in bridges.get((locus, b1, b2), {}).values()],
                        "size": abs(positions_by_key[(locus, b2)] - positions_by_key[(locus, b1)]),
                    })
        else:
            for (locus, bkp), pairs in pairs_by_key.items():
                panels.append({
                    "locus": locus, "panel": str(bkp), "class": "discordant",
                    "pairs": [p for _, p in pairs.values()
                              if _is_discordant(p, library_median, library_mad)],
                    "size": None,
                })
        for panel in panels:
            _finish_panel(panel, library_median, library_mad, args.sv_type)
        summary_rows.extend(_summary_rows(bam_name, library, library_mad, panels))

        # Pooled accumulation
        if args.plot_mode in ("pooled", "both"):
            pooled_library.extend(library)
            for panel in panels:
                pooled = pooled_panels.setdefault(
                    (panel["locus"], panel["panel"]),
                    {**panel, "pairs": [], "sizes": []},
                )
                pooled["pairs"].extend(panel["pairs"])
                pooled["sizes"].append(panel["size"])
            if canonical_pairs:
                for key in canonical_pairs:
                    b1, b2 = key
                    pooled_softclip_by_pair.setdefault(key, [])
                    pooled_softclip_by_pair[key].extend(softclip_positions_by_bkp.get(b1, []))
                    pooled_softclip_by_pair[key].extend(softclip_positions_by_bkp.get(b2, []))
            else:
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
            for i, search_region in enumerate(search_regions):
                # unplaced unmapped reads are scanned once, with the first region
                all_matches.extend(
                    find_reads_with_multiple_softclips_aho(
                        bam_file, automaton, search_region, allowed_pairs,
                        include_unplaced=(i == 0),
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
        if args.plot_insert and args.plot_mode in ("per-bam", "both") and panels:
            _plot_insert_sizes(panels, f"{args.output}_{bam_name}_insert_dist.png")

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
    if args.plot_mode in ("pooled", "both") and pooled_panels:
        library_median, library_mad = _median_mad(pooled_library)
        panels = list(pooled_panels.values())
        for panel in panels:
            if panel["class"] == "bridging":
                panel["size"] = float(np.median(panel["sizes"]))
            _finish_panel(panel, library_median, library_mad, args.sv_type)
        summary_rows.extend(_summary_rows("pooled", pooled_library, library_mad, panels))
        if args.plot_insert:
            print("Plotting pooled insert-size distributions by pair...")
            _plot_insert_sizes(panels, f"{args.output}_pooled_insert_dist_by_pair.png")
        if args.plot_softclip and pooled_softclip_by_pair:
            print("Plotting pooled softclip distributions by pair...")
            _plot_softclip_positions(
                pooled_softclip_by_pair, bkp_to_expected_pos,
                f"{args.output}_pooled_softclip_dist_by_pair.png",
                bam_name="pooled",
                pair_keys=list(pooled_softclip_by_pair.keys()),
            )

    summary_columns = ["bam", "locus", "panel", "class", "orientation", "n", "median",
                       "q1", "q3", "mad", "sv_type", "expected"]
    pd.DataFrame(summary_rows, columns=summary_columns).to_csv(
        summary_output, sep="\t", index=False, na_rep="NA", float_format="%.1f"
    )
    print(f"\nBridging pairs saved to {bridges_output}")
    print(f"Insert size summary saved to {summary_output}")
    print("Done.")
