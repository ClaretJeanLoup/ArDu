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
import matplotlib.colors
import matplotlib.patches
import matplotlib.pyplot as plt
import matplotlib.ticker
import numpy as np
import pandas as pd
import pysam

matplotlib.use("Agg")

from ardu import bam as bam_utils
from ardu.bam import passes_read_filters, read_pair_info
from ardu.junction_align import (
    EXPECTED,
    align_junction,
    expected_junction,
    mafft_available,
    microhomology_split,
)
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
    parser.add_argument("--plot-evidence", action="store_true",
                        help="Plot, per sample and breakpoint pair, the softclips at both "
                             "breakpoints and the bridging pair inserts, with the evidence "
                             "summary in the title. With --junction and MAFFT installed, "
                             "also align the junction and split reads and add their "
                             "consensus around the junction.")
    parser.add_argument("--reference", metavar="FASTA",
                        help="Indexed reference FASTA, used for the expected junction "
                             "sequence and microhomology. Defaults to the --blast reference.")
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
    window: dict,
) -> None:
    """Extract the softclips of the reads in *region* into *window*.

    Also writes to the shared FASTA / TSV outputs. *window* collects
    ``seqs`` (clip sequences), ``clips`` (``(position, side)``, side ``L`` for
    a left clip at the first aligned base, ``R`` for a right clip at the last
    one), ``pairs`` (``{read_name: (read_start, PairInfo)}``) and ``splits``
    (split reads, see :func:`_split_read`).
    """
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
                window["pairs"][read.query_name] = (read.reference_start, pair)
            insert_size = pair.insert if pair is not None else None
            split = _split_read(read)
            if split is not None:
                window["splits"][(read.query_name, read.is_read1)] = split

            # LEFT softclip
            if read.cigar[0][0] == 4 and read.cigar[0][1] >= min_clip:
                left_clip = read.query_sequence[: read.cigar[0][1]]
                window["seqs"].append(left_clip)
                fasta_out.write(
                    f">{chrom}:{mapping_position}_{read.query_name}_{breakpoint}_L\n"
                    f"{left_clip}\n"
                )
                tsv_out.write(
                    f"{bam_file}\t{read.query_name}\t{mapping_position}\t"
                    f"{mapping_position}\tleftmost_clip\t{cigar}\t{mate_position}\t"
                    f"{insert_size}\t{region}\t{breakpoint}\n"
                )
                window["clips"].append((mapping_position, "L"))

            # RIGHT softclip
            if read.cigar[-1][0] == 4 and read.cigar[-1][1] >= min_clip:
                right_clip = read.query_sequence[-read.cigar[-1][1]:]
                clip_position = read.reference_end
                window["seqs"].append(right_clip)
                fasta_out.write(
                    f">{chrom}:{clip_position}_{read.query_name}_{breakpoint}_R\n"
                    f"{right_clip}\n"
                )
                tsv_out.write(
                    f"{bam_file}\t{read.query_name}\t{mapping_position}\t"
                    f"{clip_position}\trightmost_clip\t{cigar}\t{mate_position}\t"
                    f"{insert_size}\t{region}\t{breakpoint}\n"
                )
                window["clips"].append((clip_position, "R"))


# ---------------------------------------------------------------------------
# Split reads and junction geometry
# ---------------------------------------------------------------------------

_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def _clipped_end(start: int, cigar: str) -> tuple[int, str] | None:
    """``(position, side)`` of the clipped end of an alignment starting at
    1-based *start*: the first aligned base for a left clip (``L``), the last
    one for a right clip (``R``). The longer clip wins; ``None`` if unclipped.
    """
    ops = [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]
    lead = trail = 0
    for n, op in ops:
        if op not in "SH":
            break
        lead += n
    for n, op in reversed(ops):
        if op not in "SH":
            break
        trail += n
    if lead == trail == 0:
        return None
    if lead >= trail:
        return start, "L"
    return start + sum(n for n, op in ops if op in "MDN=X") - 1, "R"


def _soft_clip_length(cigartuples) -> int:
    """Soft-clipped bases at the start of *cigartuples*, past any hard clip."""
    n = 0
    for op, length in cigartuples:
        if op == 5:
            continue
        if op != 4:
            break
        n += length
    return n


def _split_read(read) -> tuple | None:
    """Junction ends of a split read: ``(position, side, [(chrom, position, side), ...],
    sequence, offset)`` for the primary alignment and each supplementary
    alignment of its SA tag; *offset* is the index of the first base right of
    the junction in the read sequence.
    """
    if (read.is_supplementary or not read.has_tag("SA")
            or not passes_read_filters(read)):
        return None
    here = _clipped_end(read.reference_start + 1, read.cigarstring)
    if here is None:
        return None
    others = []
    for entry in read.get_tag("SA").rstrip(";").split(";"):
        chrom, pos, _strand, cigar, mapq, _nm = entry.split(",")
        end = _clipped_end(int(pos), cigar)
        if end is not None and int(mapq) >= bam_utils.MIN_MAPPING_QUALITY:
            others.append((chrom, *end))
    if not others:
        return None
    if here[1] == "L":
        offset = _soft_clip_length(read.cigartuples)
    else:
        offset = read.query_length - _soft_clip_length(reversed(read.cigartuples))
    return (*here, others, read.query_sequence, offset)


def _split_junctions(windows: dict, k1: tuple, k2: tuple) -> dict[tuple, list]:
    """Split reads joining windows *k1* and *k2*, grouped by junction
    ``(pos1, side1, pos2, side2)`` with the ends at *k1* first, as
    ``{junction: [(read_id, sequence, offset), ...]}``. Each read counts once.
    """
    junctions: dict[tuple, list] = {}
    seen: set = set()
    for here, there in ((k1, k2), (k2, k1)):
        chrom, start, end = windows[there]["region"]
        for read_id, (pos, side, others, seq, offset) in windows[here]["splits"].items():
            if read_id in seen:
                continue
            for o_chrom, o_pos, o_side in others:
                if o_chrom == chrom and start < o_pos <= end:
                    seen.add(read_id)
                    junction = ((pos, side, o_pos, o_side) if here == k1
                                else (o_pos, o_side, pos, side))
                    junctions.setdefault(junction, []).append((read_id, seq, offset))
                    break
    return junctions


def _top_clip(clips: list[tuple[int, str]]) -> tuple[int, str] | None:
    """The ``(position, side)`` with the most softclips."""
    return Counter(clips).most_common(1)[0][0] if clips else None


def _sv_geometry(a: tuple[int, str], b: tuple[int, str]) -> tuple[str | None, int | None]:
    """SV type and size from the clipped ends at two breakpoints.

    A tandem duplication leaves left clips at its first base and right clips at
    its last; a deletion leaves right clips at the last base before it and left
    clips at the first base after it.
    """
    (left_pos, left_side), (right_pos, right_side) = sorted((a, b))
    if (left_side, right_side) == ("L", "R"):
        return "duplication", right_pos - left_pos + 1
    if (left_side, right_side) == ("R", "L"):
        return "deletion", right_pos - left_pos - 1
    return None, None


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
_SIDE_COLOURS = {"L": "steelblue", "R": "darkorange"}
_SIDE_LABELS = {"L": "left clips", "R": "right clips"}


def _grid(num: int):
    num = max(1, num)
    cols = min(num, 3)
    rows = math.ceil(num / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows))
    return fig, np.array(axes).reshape(-1)


def _save(output_path: str) -> None:
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_path}")


def _count_axis(ax, counts, label: str) -> None:
    # a log axis only helps when counts span orders of magnitude
    if np.max(counts) >= 20:
        ax.set_yscale("log")
        label += " (log)"
    else:
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
    ax.set_ylabel(label, fontsize=9)


def _bp_axis(ax, label: str) -> None:
    ax.set_xlabel(label, fontsize=9)
    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda x, _: f"{int(x):,}")
    )
    ax.tick_params(axis="x", labelrotation=30, labelsize=8)


def _draw_insert_panel(ax, panel: dict) -> None:
    """Insert sizes of *panel*, stacked by orientation, with the expected insert."""
    if not panel["pairs"]:
        ax.text(0.5, 0.5, "no pairs", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
    else:
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
        _count_axis(ax, counts, "Count")
        _bp_axis(ax, "Insert size (bp)")

    if panel["expected"] is not None:
        ax.axvline(
            panel["expected"], color="red", linestyle="--", linewidth=1.5,
            label=f"expected {panel['sv_type']}: {panel['expected']:,.0f} bp",
        )
    title = f"{panel['locus']} | BKP {panel['panel']} ({panel['class']} pairs)"
    if panel["library_median"] is not None:
        title += f"\nlibrary insert {panel['library_median']:,.0f} bp"
    ax.set_title(title, fontsize=10)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(fontsize=8)


def _plot_insert_sizes(panels: list[dict], output_path: str) -> None:
    """One insert panel per breakpoint pair (bridging pairs) or per breakpoint
    (discordant pairs); see :func:`_draw_insert_panel`.
    """
    _, axes = _grid(len(panels))
    for ax, panel in zip(axes, panels):
        _draw_insert_panel(ax, panel)
    for ax in axes[len(panels):]:
        ax.axis("off")
    _save(output_path)


def _draw_softclip_panel(ax, title: str, clips: list, predicted, refined) -> None:
    """Softclip positions binned over the window, stacked by clip side, with
    the breakpoint predicted from coverage and the refined junction end.
    """
    ax.set_title(title, fontsize=10)
    if not clips:
        ax.text(0.5, 0.5, "no softclips", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        return
    positions = [p for p, _ in clips]
    span = max(positions) - min(positions) + 1
    width = max(1, math.ceil(span / 100))
    bins = np.arange(min(positions), max(positions) + width + 1, width)
    sides = [side for side in _SIDE_COLOURS if any(c[1] == side for c in clips)]
    counts, _, _ = ax.hist(
        [[p for p, side_ in clips if side_ == side] for side in sides],
        bins=bins, stacked=True,
        color=[_SIDE_COLOURS[side] for side in sides],
        label=[f"{_SIDE_LABELS[side]} (n={sum(c[1] == side for c in clips)})" for side in sides],
    )
    if predicted is not None:
        ax.axvline(predicted, color="red", linewidth=1.2,
                   label=f"predicted: {predicted:,.0f}")
    if refined is not None:
        pos, side, source = refined
        offset = f", {pos - predicted:+,.0f} bp" if predicted is not None else ""
        ax.axvline(pos, color="black", linestyle="--", linewidth=1.2,
                   label=f"junction ({source}): {pos:,} {side}{offset}")
    _count_axis(ax, counts, "Softclips")
    _bp_axis(ax, f"Genomic position ({width} bp bins)")
    ax.legend(fontsize=7, loc="upper right")


def _plot_softclip_positions(windows: list[dict], output_path: str) -> None:
    """One softclip panel per breakpoint window; see :func:`_draw_softclip_panel`."""
    _, axes = _grid(len(windows))
    for ax, w in zip(axes, windows):
        _draw_softclip_panel(ax, w["title"], w["clips"], w["predicted"], w["refined"])
    for ax in axes[len(windows):]:
        ax.axis("off")
    _save(output_path)


def _evidence_text(ev: dict) -> str:
    if ev["sv_type"] and ev["size"]:
        call = (f"{ev['sv_type']} {ev['bkp1_refined']:,}-{ev['bkp2_refined']:,} "
                f"({ev['size']:,} bp)" if ev["bkp1_refined"] else
                f"{ev['sv_type']} (~{ev['size']:,} bp)")
    else:
        call = "no SV call"
    junction_reads = "NA" if ev["junction_reads"] is None else ev["junction_reads"]
    pairs = f"{ev['bridging_pairs']} bridging pairs"
    if ev["bridging_pairs"]:
        pairs += f" ({ev['pair_orientation']}, median insert {ev['insert_median']:,.0f} bp)"
    text = (f"{ev['bam']} | {ev['locus']} BKP {ev['pair']}: {call}\n"
            f"{ev['split_reads']} split reads, {junction_reads} junction reads, {pairs}")
    if ev.get("aligned_reads"):
        text += f"\n{ev['aligned_reads']} reads aligned"
        if ev.get("junction_agreement") is not None:
            text += f", {ev['junction_agreement']:.0%} agreement at the junction"
        if ev.get("inserted_bases") is not None:
            text += f", inserted bases: {ev['inserted_bases'] or 'none'}"
    if ev.get("microhomology") is not None:
        text += f", microhomology: {ev['microhomology'] or 'none'}"
    return text


_BASE_COLOURS = {"A": "#2ca02c", "C": "#1f77b4", "G": "#ff7f0e", "T": "#d62728"}


def _draw_alignment_panel(ax, alignment: dict, flank: int = 30, max_reads: int = 40) -> None:
    """Aligned reads around the junction: grey where a read agrees with the
    consensus, coloured by base where it does not, dark for gaps.
    """
    aligned, columns, junction_col = (
        alignment["aligned"], alignment["columns"], alignment["junction_col"]
    )
    covered = [c for c, (_, depth, _) in enumerate(columns) if depth]
    centre = junction_col if junction_col is not None else covered[len(covered) // 2]
    lo, hi = max(covered[0], centre - flank), min(covered[-1] + 1, centre + flank)

    def first_base(row):
        return len(row) - len(row.lstrip("-"))

    reads = sorted((n for n in aligned if n != EXPECTED), key=lambda n: first_base(aligned[n]))
    names = ([EXPECTED] if EXPECTED in aligned else []) + reads[:max_reads]
    codes = np.zeros((len(names), hi - lo))
    for i, name in enumerate(names):
        row = aligned[name]
        start, end = first_base(row), len(row.rstrip("-")) - 1
        for j, c in enumerate(range(lo, hi)):
            if name != EXPECTED and not start <= c <= end:
                continue
            if row[c] == columns[c][0]:
                codes[i, j] = 1
            elif row[c] == "-":
                codes[i, j] = 2
            else:
                codes[i, j] = 3 + "ACGT".find(row[c]) if row[c] in "ACGT" else 2
    cmap = matplotlib.colors.ListedColormap(
        ["white", "#d9d9d9", "#555555", *_BASE_COLOURS.values()]
    )
    ax.imshow(codes, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=6)
    ax.set_xticks(range(hi - lo))
    ax.set_xticklabels([columns[c][0] or " " for c in range(lo, hi)], fontsize=6,
                       family="monospace")
    ax.set_xlabel("Consensus", fontsize=9)
    ax.set_yticks([0] if names and names[0] == EXPECTED else [])
    ax.set_yticklabels(["expected"] if names and names[0] == EXPECTED else [], fontsize=8)
    if junction_col is not None:
        ax.axvline(junction_col - lo - 0.5, color="black", linewidth=1.5)
    stats = alignment["stats"]
    title = f"{stats['aligned_reads']} reads aligned (MAFFT)"
    if stats["junction_depth"] is not None:
        title += (f", junction depth {stats['junction_depth']}, "
                  f"agreement {stats['junction_agreement']:.0%}")
    if len(reads) > max_reads:
        title += f"\nfirst {max_reads} reads shown"
    ax.set_title(title, fontsize=10)
    ax.legend(
        handles=[matplotlib.patches.Patch(color=c, label=b) for b, c in _BASE_COLOURS.items()]
        + [matplotlib.patches.Patch(color="#555555", label="gap")],
        fontsize=7, loc="upper right", ncol=5,
    )


_SEQUENCE_COLOURS = {"left": "steelblue", "right": "darkorange", "inserted": "red",
                     "microhomology": "mediumpurple", "low support": "#999999"}


def _draw_sequence_panel(ax, alignment: dict, ev: dict, line_length: int = 100) -> None:
    """The consensus sequence, wrapped, coloured by side of the junction, with
    inserted bases, microhomology and weakly supported (lower case) bases marked.
    """
    seq, junction = alignment["consensus"], alignment["consensus_junction"]
    inserted = len(alignment["stats"]["inserted_bases"] or "")
    mh_left, mh_right = (len(x) for x in ev.get("_microhomology") or ("", ""))
    n_lines = math.ceil(len(seq) / line_length)

    def category(i: int) -> str:
        if junction is None:
            return "left"
        if junction - inserted <= i < junction and inserted:
            return "inserted"
        if junction - inserted - mh_left <= i < junction - inserted or junction <= i < junction + mh_right:
            return "microhomology"
        return "left" if i < junction else "right"

    for i, base in enumerate(seq):
        x, y = i % line_length, -(i // line_length)
        colour = _SEQUENCE_COLOURS["low support"] if base.islower() else _SEQUENCE_COLOURS[category(i)]
        ax.text(x, y, base.upper(), family="monospace", fontsize=9, ha="center", va="center",
                color=colour, fontweight="bold" if category(i) != "left" and category(i) != "right" else None)
    for line in range(n_lines):
        ax.text(-1.5, -line, f"{line * line_length + 1}", ha="right", va="center", fontsize=7,
                color="#555555")
    if junction is not None:
        x, y = junction % line_length, -(junction // line_length)
        ax.plot([x - 0.5, x - 0.5], [y - 0.45, y + 0.45], color="black", linewidth=1.5)

    caption = ""
    junction_ends = {side: pos for pos, side in ev["_ends"]} if ev.get("_ends") else {}
    if junction is not None and set(junction_ends) == {"L", "R"}:
        left_len = junction - inserted
        right_len = len(seq) - junction
        caption = (f"left of |: {ev['chrom']}:{junction_ends['R'] - left_len + 1:,}-"
                   f"{junction_ends['R']:,}    right of |: {ev['chrom']}:"
                   f"{junction_ends['L']:,}-{junction_ends['L'] + right_len - 1:,}")
    ax.set_xlim(-6, line_length + 1)
    ax.set_ylim(-n_lines + 0.4, 0.8)
    ax.axis("off")
    ax.set_title(f"Junction consensus ({len(seq)} bp)" + (f"\n{caption}" if caption else ""),
                 fontsize=10)
    shown = {"left", "right", "low support"}
    if inserted:
        shown.add("inserted")
    if mh_left or mh_right:
        shown.add("microhomology")
    ax.legend(
        handles=[matplotlib.patches.Patch(color=c, label=label)
                 for label, c in _SEQUENCE_COLOURS.items() if label in shown],
        fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 0.0), ncol=len(shown),
        frameon=False,
    )


def _plot_evidence(ev: dict, windows: tuple[dict, dict], panel: dict, output_path: str,
                   alignment: dict | None = None) -> None:
    """Softclips at both breakpoints, bridging pair inserts and, when the
    junction reads were aligned, the alignment around the junction and the
    consensus sequence, with the evidence summarised in the title.
    """
    if alignment is None:
        fig, axes = plt.subplots(1, 3, figsize=(18, 4.5))
    else:
        n_lines = math.ceil(len(alignment["consensus"]) / 100)
        fig = plt.figure(figsize=(14, 9 + 0.35 * n_lines + 1))
        grid = fig.add_gridspec(3, 2, height_ratios=[4.5, 4.5, 0.35 * n_lines + 1])
        axes = [fig.add_subplot(grid[r, c]) for r in (0, 1) for c in (0, 1)]
        axes.append(fig.add_subplot(grid[2, :]))
    for ax, w in zip(axes, windows):
        _draw_softclip_panel(ax, w["title"], w["clips"], w["predicted"], w["refined"])
    _draw_insert_panel(axes[2], panel)
    if alignment is not None:
        _draw_alignment_panel(axes[3], alignment)
        _draw_sequence_panel(axes[4], alignment, ev)
    fig.suptitle(_evidence_text(ev), fontsize=11, y=1.02)
    _save(output_path)


# ---------------------------------------------------------------------------
# Evidence per breakpoint pair
# ---------------------------------------------------------------------------

def _locus_pairs(windows: dict, canonical_pairs: list) -> list[tuple]:
    """``(locus, b1, b2)`` for the --pairs present in each locus, or every pair
    of breakpoints within a locus without --pairs.
    """
    out = []
    for locus in dict.fromkeys(k[0] for k in windows):
        bkps = sorted(b for l, b in windows if l == locus)
        candidates = canonical_pairs or [
            (a, b) for i, a in enumerate(bkps) for b in bkps[i + 1:]
        ]
        out += [(locus, b1, b2) for b1, b2 in candidates if b1 in bkps and b2 in bkps]
    return out


def _junction_search(bam_file: str, windows: dict, locus_pairs: list) -> tuple[list, dict]:
    """Reads containing softclips from both breakpoints of a pair, and their
    number per ``(locus, b1, b2)``.
    """
    allowed = set()
    for locus, b1, b2 in locus_pairs:
        allowed |= {((locus, b1), (locus, b2)), ((locus, b2), (locus, b1))}
    keys = sorted({k for pair in allowed for k in pair})
    if not any(windows[k]["seqs"] for k in keys):
        return [], {}
    automaton = build_softclip_automaton({k: windows[k]["seqs"] for k in keys})

    # one fetch over the span of each locus when its windows share a chromosome
    search_regions = []
    for locus in dict.fromkeys(k[0] for k in keys):
        regions = [windows[k]["region"] for k in keys if k[0] == locus]
        if len({r[0] for r in regions}) == 1:
            search_regions.append((regions[0][0], min(r[1] for r in regions),
                                   max(r[2] for r in regions)))
        else:
            search_regions.extend(regions)

    matches, seen = [], set()
    for i, region in enumerate(search_regions):
        # unplaced unmapped reads are scanned once, with the first region
        for m in find_reads_with_multiple_softclips_aho(
            bam_file, automaton, region, allowed, include_unplaced=(i == 0),
        ):
            key = (m["read_name"], m["breakpoints"])
            if key not in seen:
                seen.add(key)
                matches.append(m)

    names: dict[tuple, set] = {}
    for m in matches:
        (locus, a), (_, b) = m["breakpoints"]
        names.setdefault((locus, *_canonical_pair_key(a, b)), set()).add(m["read_name"])
    return matches, {k: len(v) for k, v in names.items()}


def _pair_evidence(
    bam_name: str, locus, b1: int, b2: int, windows: dict, bridges: dict,
    junction_counts: dict | None, library_median, library_mad, sv_type: str,
) -> tuple[dict, dict]:
    """Evidence row and bridging insert panel for breakpoints *b1*, *b2* of *locus*.

    The junction is the one most split reads support, or the top softclip of
    each window without split reads. The SV type is, in order: --sv-type, the
    split read geometry, the bridging pair orientation, the softclip geometry.
    """
    k1, k2 = (locus, b1), (locus, b2)
    w1, w2 = windows[k1], windows[k2]
    junctions = _split_junctions(windows, k1, k2)
    ends, refined_by, split_reads = (None, None), None, []
    if junctions:
        (p1, s1, p2, s2), split_reads = max(junctions.items(), key=lambda j: len(j[1]))
        ends, refined_by = ((p1, s1), (p2, s2)), "split_reads"
    elif w1["clips"] and w2["clips"]:
        ends, refined_by = (_top_clip(w1["clips"]), _top_clip(w2["clips"])), "softclips"
    geometry_type, size = _sv_geometry(*ends) if refined_by else (None, None)

    chosen = sv_type if sv_type != "auto" else (
        geometry_type if refined_by == "split_reads" else None
    )
    pairs = [p for _, p in bridges.get((locus, b1, b2), {}).values()]
    panel = {
        "locus": locus, "panel": f"{b1}-{b2}", "class": "bridging", "pairs": pairs,
        "size": size if size is not None else abs(w2["predicted"] - w1["predicted"]),
    }
    _finish_panel(panel, library_median, library_mad, chosen or "auto")

    ev = {
        "bam": bam_name, "locus": locus, "pair": f"{b1}-{b2}", "chrom": w1["region"][0],
        "bkp1_predicted": w1["predicted"], "bkp2_predicted": w2["predicted"],
        "bkp1_refined": ends[0][0] if ends[0] else None,
        "bkp1_side": ends[0][1] if ends[0] else None,
        "bkp2_refined": ends[1][0] if ends[1] else None,
        "bkp2_side": ends[1][1] if ends[1] else None,
        "refined_by": refined_by,
        "sv_type": chosen or panel["sv_type"] or geometry_type,
        "size": size,
        "split_reads": sum(len(v) for v in junctions.values()),
        "split_reads_at_junction": len(split_reads),
        "junction_reads": (None if junction_counts is None
                           else junction_counts.get((locus, b1, b2), 0)),
        "softclips_bkp1": len(w1["clips"]),
        "softclips_bkp2": len(w2["clips"]),
        "bridging_pairs": len(pairs),
        "pair_orientation": (Counter(p.orientation for p in pairs).most_common(1)[0][0]
                             if pairs else None),
        "insert_median": float(np.median([p.insert for p in pairs])) if pairs else None,
        "insert_expected": panel["expected"],
        # kept for the junction alignment, not written to the summary
        "_ends": ends if refined_by else None,
        "_split_reads": split_reads,
    }
    return ev, panel


def _junction_alignment(ev: dict, matches: list, k1: tuple, k2: tuple, fasta,
                        prefix: str) -> dict | None:
    """Align the junction reads and the split reads of one breakpoint pair with
    MAFFT, write the alignment, and add its statistics to *ev*.
    """
    reads: dict[str, str] = {}
    split_offsets: dict[str, int] = {}
    for (name, is_read1), seq, offset in ev["_split_reads"]:
        read_id = f"{name}/{1 if is_read1 else 2}"
        reads[read_id] = seq
        split_offsets[read_id] = offset
    known = set(reads.values())
    for m in matches:
        if set(m["breakpoints"]) == {k1, k2} and m["sequence"] not in known:
            name = m["read_name"]
            while name in reads:
                name += "'"
            reads[name] = m["sequence"]
            known.add(m["sequence"])
    if len(reads) < 2:
        return None

    expected = (expected_junction(fasta, ev["chrom"], ev["_ends"])
                if fasta is not None and ev["_ends"] else None)
    alignment = align_junction(reads, split_offsets, expected)
    ev.update(alignment["stats"])

    with open(f"{prefix}_junction_alignment.txt", "w") as fh:
        fh.writelines(f"# {line}\n" for line in _evidence_text(ev).splitlines())
        fh.write(f"# consensus: {alignment['consensus']}\n")
        fh.write(alignment["text"])
    with open(f"{prefix}_junction_msa.fasta", "w") as fh:
        consensus_row = "".join(symbol or "-" for symbol, _, _ in alignment["columns"])
        fh.write(f">consensus\n{consensus_row}\n")
        for name, row in alignment["aligned"].items():
            fh.write(f">{name}\n{row}\n")
    print(f"Saved {prefix}_junction_alignment.txt")
    return alignment


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
    evidence_output = f"{args.output}_junction_summary.tsv"

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

    with open(tsv_output, "w") as out:
        out.write(
            "bam_file\tRead_name\tPosition\tClipping_position\t"
            "Clipping_end\tCIGAR\tMate_position\tInsert_size\tRegion\tBreakpoint\n"
        )
    with open(bridges_output, "w") as bf:
        bf.write("bam_file\tlocus\tread_name\tread_pos\tmate_pos\tinsert_size\t"
                 "orientation\tproper_pair\tbreakpoint_bridge\n")

    summary_rows: list[dict] = []
    evidence_rows: list[dict] = []
    pooled_library: list[int] = []
    pooled_panels: dict = {}
    pooled_windows: dict = {}
    output_dir = os.path.dirname(args.output) or "."
    reference = getattr(args, "reference", None) or args.blast
    fasta = pysam.FastaFile(reference) if reference else None
    align_reads = args.junction and args.plot_evidence
    if align_reads and not mafft_available():
        print("Warning: mafft not found; junction reads will not be aligned.")
        align_reads = False

    # ------------------------------------------------------------------
    # Per-BAM processing
    # ------------------------------------------------------------------
    for bam_basename, bam_file in bam_mapping.items():
        print(f"\n=== PROCESSING BAM: {bam_basename} ===")
        bam_name = bam_basename.removesuffix(".bam")

        # Breakpoint numbers restart for each sample and locus, so windows are
        # keyed by (locus, breakpoint) within this BAM.
        windows: dict[tuple, dict] = {}
        bkp_subset = breakpoints_df[breakpoints_df["bam_file"] == bam_basename]
        for _, row in bkp_subset.iterrows():
            bkp = row["breakpoint_number"]
            pos = row["position"]
            key = (row.get("locus", ""), bkp)
            windows[key] = {
                "region": (row["chromosome"], pos - args.extension, pos + args.extension),
                "predicted": pos, "seqs": [], "clips": [], "pairs": {}, "splits": {},
            }
            _extract_softclips(bam_file, windows[key]["region"], args.size,
                               fasta_output, tsv_output, bkp, windows[key])

        # ---- Insert sizes ----
        pairs_by_key = {k: w["pairs"] for k, w in windows.items()}
        library = _library_inserts(pairs_by_key)
        library_median, library_mad = _median_mad(library)
        if library_median is None:
            print(f"No proper FR pairs in the windows of {bam_name}; "
                  "discordant pairs are those not flagged proper by the aligner.")
        elif 2 * args.extension < library_median:
            print(f"Warning: windows of {2 * args.extension} bp are shorter than the "
                  f"library insert ({library_median:,.0f} bp) of {bam_name}; few pairs "
                  "will bridge two breakpoints. Consider a larger -e.")

        bridges = _bridging_pairs(
            pairs_by_key, {k: w["region"] for k, w in windows.items()},
            library_median, library_mad,
        )
        with open(bridges_output, "a") as bf:
            for (locus, b1, b2), pairs in bridges.items():
                chrom = windows[(locus, b1)]["region"][0]
                for name, (read_start, pair) in pairs.items():
                    bf.write(
                        f"{bam_basename}\t{locus}\t{name}\t{chrom}:{read_start + 1}\t"
                        f"{chrom}:{pair.mate_start + 1}\t{pair.insert}\t"
                        f"{pair.orientation}\t{pair.proper}\t{b1}-{b2}\n"
                    )

        # ---- Junction search ----
        locus_pairs = _locus_pairs(windows, canonical_pairs)
        junction_counts = None
        matches: list = []
        if args.junction:
            print(f"Searching junction reads for {bam_name}...")
            matches, junction_counts = _junction_search(bam_file, windows, locus_pairs)
            if getattr(args, "junction_fastq", False):
                write_junction_reads_fastq(bam_file, matches, output_dir, bam_name)
            else:
                write_junction_reads_fasta(matches, output_dir, bam_name)

        # ---- Evidence per breakpoint pair ----
        refined_by_key: dict[tuple, tuple] = {}
        pair_results = []
        for locus, b1, b2 in locus_pairs:
            ev, panel = _pair_evidence(
                bam_name, locus, b1, b2, windows, bridges, junction_counts,
                library_median, library_mad, args.sv_type,
            )
            alignment = None
            if fasta is not None and ev["_ends"]:
                ev["_microhomology"] = microhomology_split(fasta, ev["chrom"], ev["_ends"])
                if ev["_microhomology"] is not None:
                    ev["microhomology"] = "".join(ev["_microhomology"])
            if align_reads:
                alignment = _junction_alignment(
                    ev, matches, (locus, b1), (locus, b2), fasta,
                    f"{args.output}_{bam_name}_{locus}_{b1}-{b2}",
                )
            evidence_rows.append(ev)
            pair_results.append((locus, b1, b2, ev, panel, alignment))
            for bkp, n in ((b1, 1), (b2, 2)):
                if ev[f"bkp{n}_refined"] is not None:
                    refined_by_key.setdefault((locus, bkp), (
                        ev[f"bkp{n}_refined"], ev[f"bkp{n}_side"], ev["refined_by"]
                    ))
        for key, w in windows.items():
            if key not in refined_by_key and w["clips"]:
                refined_by_key[key] = (*_top_clip(w["clips"]), "softclips")
        window_plots = {
            key: {"title": f"{bam_name} | {key[0]} BKP {key[1]}", "clips": w["clips"],
                  "predicted": w["predicted"], "refined": refined_by_key.get(key)}
            for key, w in windows.items()
        }

        if canonical_pairs:
            panels = [panel for *_, panel, _ in pair_results]
        else:
            panels = []
            for (locus, bkp), pairs in pairs_by_key.items():
                panel = {
                    "locus": locus, "panel": str(bkp), "class": "discordant",
                    "pairs": [p for _, p in pairs.values()
                              if _is_discordant(p, library_median, library_mad)],
                    "size": None,
                }
                _finish_panel(panel, library_median, library_mad, args.sv_type)
                panels.append(panel)
        summary_rows.extend(_summary_rows(bam_name, library, library_mad, panels))

        # ---- Per-BAM plots ----
        if args.plot_insert and args.plot_mode in ("per-bam", "both") and panels:
            _plot_insert_sizes(panels, f"{args.output}_{bam_name}_insert_dist.png")
        if args.plot_softclip and args.plot_mode in ("per-bam", "both") and windows:
            _plot_softclip_positions(list(window_plots.values()),
                                     f"{args.output}_{bam_name}_softclip_dist.png")
        if args.plot_evidence:
            for locus, b1, b2, ev, panel, alignment in pair_results:
                _plot_evidence(
                    ev, (window_plots[(locus, b1)], window_plots[(locus, b2)]), panel,
                    f"{args.output}_{bam_name}_{locus}_{b1}-{b2}_evidence.png",
                    alignment,
                )

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
            for key, w in window_plots.items():
                pooled = pooled_windows.setdefault(
                    key, {"clips": [], "predicted": [], "refined": []}
                )
                pooled["clips"].extend(w["clips"])
                pooled["predicted"].append(w["predicted"])
                if w["refined"] is not None:
                    pooled["refined"].append(w["refined"])

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
    if args.plot_mode in ("pooled", "both") and args.plot_softclip and pooled_windows:
        print("Plotting pooled softclip distributions...")
        _plot_softclip_positions(
            [{"title": f"pooled | {locus} BKP {bkp}", "clips": w["clips"],
              "predicted": float(np.median(w["predicted"])),
              "refined": Counter(w["refined"]).most_common(1)[0][0] if w["refined"] else None}
             for (locus, bkp), w in pooled_windows.items()],
            f"{args.output}_pooled_softclip_dist.png",
        )

    summary_columns = ["bam", "locus", "panel", "class", "orientation", "n", "median",
                       "q1", "q3", "mad", "sv_type", "expected"]
    pd.DataFrame(summary_rows, columns=summary_columns).to_csv(
        summary_output, sep="\t", index=False, na_rep="NA", float_format="%.1f"
    )
    evidence_columns = [
        "bam", "locus", "pair", "chrom", "bkp1_predicted", "bkp2_predicted",
        "bkp1_refined", "bkp1_side", "bkp2_refined", "bkp2_side", "refined_by",
        "sv_type", "size", "split_reads", "split_reads_at_junction", "junction_reads",
        "softclips_bkp1", "softclips_bkp2", "bridging_pairs", "pair_orientation",
        "insert_median", "insert_expected", "microhomology", "aligned_reads",
        "consensus_length", "junction_depth", "junction_agreement", "inserted_bases",
        "expected_identity",
    ]
    pd.DataFrame(evidence_rows, columns=evidence_columns).astype(
        {c: "Int64" for c in ("bkp1_refined", "bkp2_refined", "size", "junction_reads",
                              "aligned_reads", "consensus_length", "junction_depth")}
    ).to_csv(evidence_output, sep="\t", index=False, na_rep="NA", float_format="%.0f")
    print(f"\nBridging pairs saved to {bridges_output}")
    print(f"Insert size summary saved to {summary_output}")
    print(f"Junction evidence saved to {evidence_output}")
    print("Done.")
