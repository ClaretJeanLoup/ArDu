"""
``ardu softclips`` — extract soft-clipped sequences and optionally run BLAST.

Corresponds to ArDu_softclips.py, refactored as a subcommand with a
``register()`` / ``run()`` interface.
"""

from __future__ import annotations

import math
import os
import re
import subprocess
from collections import Counter

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pysam
import seaborn as sns

matplotlib.use("Agg")


# ---------------------------------------------------------------------------
# Argparse registration
# ---------------------------------------------------------------------------

def register(subparsers) -> None:
    """Add the ``softclips`` sub-command to *subparsers*."""
    parser = subparsers.add_parser(
        "softclips",
        help="Extract soft-clipped sequences around breakpoints and run BLAST.",
        description=(
            "ArDu softclips — extract soft-clipped reads around each breakpoint, "
            "write FASTA, optionally BLAST against a reference, and filter hits."
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
                        help="Reference FASTA for BLAST; triggers BLAST run + filtering.")
    parser.add_argument("--pairs",
                        help="File listing allowed breakpoint pairs (TSV with source_bkp / target_bkp).")
    parser.add_argument("--plot", action="store_true",
                        help="Plot softclip density per sample per breakpoint.")

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Helpers (previously module-level in ArDu_softclips.py)
# ---------------------------------------------------------------------------

def _extract_breakpoint_from_qseqid(qseqid: str) -> int | None:
    match = re.search(r"_(\d+)_", qseqid)
    return int(match.group(1)) if match else None


def _filter_blast_hits_by_breakpoint(
    blast_output: str,
    breakpoints_df: pd.DataFrame,
    extension: int,
    filtered_output: str,
    pairs_file: str | None = None,
) -> None:
    """Filter BLAST hits so that each query softclip maps near a different breakpoint."""
    allowed_pairs: set[tuple[int, int]] = set()
    if pairs_file:
        pairs_df = pd.read_csv(pairs_file, sep="\t")
        for _, row in pairs_df.iterrows():
            source = int(row["source_bkp"])
            target = int(row["target_bkp"])
            allowed_pairs.add((source, target))
            allowed_pairs.add((target, source))

    regions_by_bkp: dict[int, tuple] = {}
    for _, row in breakpoints_df.iterrows():
        chrom = row["chromosome"]
        pos = row["position"]
        bkp = row["breakpoint_number"]
        regions_by_bkp[bkp] = (chrom, pos - extension, pos + extension)

    columns = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    ]
    blast_df = pd.read_csv(blast_output, sep="\t", names=columns)

    filtered_hits = []
    for _, row in blast_df.iterrows():
        qseqid = row["qseqid"]
        sseqid = str(row["sseqid"])
        sstart = min(row["sstart"], row["send"])
        send_val = max(row["sstart"], row["send"])

        source_bkp = _extract_breakpoint_from_qseqid(qseqid)
        if source_bkp is None:
            continue

        for target_bkp, (chrom, start, end) in regions_by_bkp.items():
            if target_bkp == source_bkp:
                continue
            if allowed_pairs and (source_bkp, target_bkp) not in allowed_pairs:
                continue
            if sseqid == chrom and not (send_val < start or sstart > end):
                filtered_hits.append(row)
                break

    pd.DataFrame(filtered_hits).to_csv(filtered_output, sep="\t", index=False, header=False)
    print(f"Filtered BLAST results saved to {filtered_output}")


def _run_blast(query_fasta: str, reference: str, blast_output: str) -> None:
    if not os.path.exists(f"{reference}.nhr"):
        print("Creating BLAST database...")
        subprocess.run(
            ["makeblastdb", "-in", reference, "-dbtype", "nucl"], check=True
        )
    print("Running BLAST...")
    subprocess.run(
        ["blastn", "-query", query_fasta, "-db", reference,
         "-out", blast_output, "-outfmt", "6"],
        check=True,
    )
    print(f"BLAST results saved in {blast_output}")


def _extract_soft_clipped_positions(
    bam_file: str,
    region: tuple,
    min_clip: int,
    output_fasta: str,
    output_tsv: str,
    breakpoint: int,
    clipping_positions: dict,
    breakpoint_center: int,
) -> None:
    """Extract soft-clipped sequences and append to FASTA / TSV outputs."""
    with pysam.AlignmentFile(bam_file, "rb") as bam, \
         open(output_fasta, "a") as fasta_out, \
         open(output_tsv, "a") as tsv_out:

        for read in bam.fetch(region[0], region[1], region[2]):
            if read.cigar is None:
                continue

            mapping_position = read.reference_start + 1
            cigar = read.cigarstring
            chrom = read.reference_name
            mate_position = (
                read.next_reference_start if not read.mate_is_unmapped else None
            )
            insert_size = (
                abs(mapping_position - mate_position)
                if mate_position is not None
                else None
            )

            # LEFT soft-clip
            if read.cigar[0][0] == 4 and read.cigar[0][1] >= min_clip:
                left_clip = read.query_sequence[: read.cigar[0][1]]
                fasta_out.write(
                    f">{chrom}:{mapping_position}_{read.query_name}_{breakpoint}_L\n"
                    f"{left_clip}\n"
                )
                tsv_out.write(
                    f"{bam_file}\t{read.query_name}\t{mapping_position}\t"
                    f"{mapping_position}\tleftmost_clip\t{cigar}\t{mate_position}\t"
                    f"{insert_size}\t{region}\t{breakpoint}\n"
                )
                relative_pos = mapping_position - breakpoint_center
                clipping_positions.setdefault(bam_file, {}).setdefault(
                    breakpoint, []
                ).append(relative_pos)

            # RIGHT soft-clip
            if read.cigar[-1][0] == 4 and read.cigar[-1][1] >= min_clip:
                right_clip = read.query_sequence[-read.cigar[-1][1] :]
                clip_position = (
                    mapping_position
                    + (read.query_length - read.cigar[-1][1])
                )
                fasta_out.write(
                    f">{chrom}:{clip_position}_{read.query_name}_{breakpoint}_R\n"
                    f"{right_clip}\n"
                )
                tsv_out.write(
                    f"{bam_file}\t{read.query_name}\t{mapping_position}\t"
                    f"{clip_position}\trightmost_clip\t{cigar}\t{mate_position}\t"
                    f"{insert_size}\t{region}\t{breakpoint}\n"
                )
                relative_pos = clip_position - breakpoint_center
                clipping_positions.setdefault(bam_file, {}).setdefault(
                    breakpoint, []
                ).append(relative_pos)


def _plot_softclip_density(clipping_positions: dict, output_prefix: str) -> None:
    samples = list(clipping_positions.keys())
    breakpoints = sorted(
        {bp for sample in clipping_positions.values() for bp in sample}
    )
    n_rows, n_cols = len(samples), len(breakpoints)
    plt.figure(figsize=(5 * n_cols, 4 * n_rows))
    plot_idx = 1

    for sample in samples:
        for bp in breakpoints:
            plt.subplot(n_rows, n_cols, plot_idx)
            positions = clipping_positions[sample].get(bp, [])
            if positions:
                sns.kdeplot(positions, fill=True)
                counts = Counter(positions)
                mode_pos = counts.most_common(1)[0][0]
                plt.axvline(mode_pos, linestyle="--")
                plt.text(
                    mode_pos,
                    plt.ylim()[1] * 0.8,
                    f"mode={mode_pos}",
                    rotation=90,
                    verticalalignment="center",
                )
            plt.title(f"{os.path.basename(sample)} | BP {bp}")
            plt.xlabel("Position relative to breakpoint")
            plt.ylabel("Density")
            plot_idx += 1

    plt.tight_layout()
    plt.savefig(f"{output_prefix}_softclip_density.png")
    plt.close()


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run(args) -> None:
    """Entry point for ``ardu softclips``."""
    fasta_output = f"{args.output}_soft_clipped_sequences.fasta"
    blast_output = f"{args.output}_BLAST_results.txt"
    tsv_output = f"{args.output}.tsv"

    breakpoints_df = pd.read_csv(args.input, sep="\t")
    clipping_positions: dict = {}

    with open(tsv_output, "w") as out:
        out.write(
            "bam_file\tRead_name\tPosition\tClipping_position\t"
            "Clipping_end\tCIGAR\tMate_position\tInsert_size\tRegion\tBreakpoint\n"
        )

    with open(args.bam_list) as fh:
        bam_files = [l.strip() for l in fh if l.strip()]

    for _, row in breakpoints_df.iterrows():
        bam_basename = row["bam_file"]
        # Resolve full path: match against bam_files list
        bam_file = next(
            (b for b in bam_files if os.path.basename(b).removesuffix(".bam") == bam_basename),
            bam_basename + ".bam",
        )
        chromosome = row["chromosome"]
        position = row["position"]
        breakpoint_number = row["breakpoint_number"]
        region = (chromosome, position - args.extension, position + args.extension)

        _extract_soft_clipped_positions(
            bam_file, region, args.size,
            fasta_output, tsv_output,
            breakpoint_number, clipping_positions, position,
        )

    if args.blast:
        _run_blast(fasta_output, args.blast, blast_output)
        filtered_output = f"{args.output}_BLAST_filtered_pairs.txt"
        _filter_blast_hits_by_breakpoint(
            blast_output, breakpoints_df, args.extension,
            filtered_output, getattr(args, "pairs", None),
        )

    if args.plot and clipping_positions:
        _plot_softclip_density(clipping_positions, args.output)
