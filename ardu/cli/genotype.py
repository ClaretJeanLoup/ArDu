"""
``ardu genotype`` — amino acids at known codons, and their copies in duplicated loci.
"""

from __future__ import annotations

import math
import os
import sys
import time

import pandas as pd
import pysam

from ardu import bam
from ardu.codons import (
    CLASSES,
    aa_class,
    allocate_copies,
    genotype_string,
    read_codons,
    read_mutations,
    reference_codon,
    total_copies,
    translate,
)


# ---------------------------------------------------------------------------
# Argparse registration
# ---------------------------------------------------------------------------

def register(subparsers) -> None:
    """Add the ``genotype`` sub-command to *subparsers*."""
    parser = subparsers.add_parser(
        "genotype",
        help="Amino acids read at known codons, and their copies given the RCN.",
        description=(
            "ArDu genotype — count the codons read over each mutation, translate "
            "them and, with the coverage output of ardu coverage, share the copies "
            "of the locus between the susceptible (S), resistant (R) and other (O) "
            "amino acids, e.g. 2S1R."
        ),
    )

    parser.add_argument("-b", "--bam", required=True,
                        help="File listing BAM paths, one per line.")
    parser.add_argument("-m", "--mutation", required=True,
                        help="Tab-delimited file: chromosome, start, stop (1-based, inclusive, "
                             "the 3 bases of the codon), gene, mutation name (G119S, L1014S/F), "
                             "strand (forward/reverse or +/-). The gene must match a locus of "
                             "the coverage file.")
    parser.add_argument("-f", "--fasta", required=True,
                        help="Reference fasta (indexed) of the bams. The reference codon must "
                             "translate to the first amino acid of each mutation name.")
    parser.add_argument("-o", "--outfile", required=True, metavar="PREFIX",
                        help="Prefix of the output file <prefix>_codons.tsv.")
    parser.add_argument("-c", "--coverage",
                        help="Coverage output of ardu coverage (<prefix>_coverage.tsv) with the "
                             "RCN of each gene. Without it every locus is taken as not duplicated.")

    parser.add_argument("--translation-table", type=int, default=1,
                        help="NCBI translation table. Default=1 (standard code).")
    parser.add_argument("--min-mapq", type=int, default=40,
                        help="Minimum mapping quality of the reads. Default=40.")
    parser.add_argument("--min-baseq", type=int, default=20,
                        help="Minimum base quality of each base of the codon. Default=20.")
    parser.add_argument("--min-depth", type=int, default=10,
                        help="Minimum number of codons read; below it the genotype is NA. "
                             "Default=10.")
    parser.add_argument("--min-reads", type=int, default=3,
                        help="Minimum reads for an amino acid to get a copy. Default=3.")
    parser.add_argument("--min-frac", type=float, default=0.05,
                        help="Minimum share of the reads for an amino acid to get a copy. "
                             "The share of half a copy, 0.5 / copies, applies when higher. "
                             "Default=0.05.")
    parser.add_argument("--dup-rcn", type=float, default=1.4,
                        help="Copies per expected copy above which a locus counts as duplicated "
                             "(RCN for a diploid locus normalised on a diploid one). Below it the "
                             "locus keeps its ploidy. Default=1.4.")
    parser.add_argument("--haploid", nargs="+", metavar="CHROM", default=[],
                        help="Chromosomes with one copy (e.g. X in males), or 'all' for a haploid "
                             "genome. The normalisation locus is taken as diploid unless 'all'.")

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _read_rcn(path: str) -> dict[str, dict[str, float | None]]:
    """``{locus: {sample: RCN}}`` from an ardu coverage output; ``None`` where RCN is NA."""
    rcn: dict[str, dict[str, float | None]] = {}
    with open(path) as fh:
        samples = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if samples is None:
                samples = parts[1:]
                continue
            rcn[parts[0]] = {}
            for sample, cell in zip(samples, parts[1:]):
                fields = cell.split(";")
                try:
                    value = float(fields[4])
                except (IndexError, ValueError):
                    value = math.nan
                rcn[parts[0]][sample] = None if math.isnan(value) else value
    return rcn


def _genotype_cell(mutation, counts, n_indel, rcn, args) -> str:
    """``Genotype;ReadsS;ReadsR;ReadsO;Depth;RCN;Codons;Flags`` for one sample."""
    reads_aa: dict[str, int] = {}
    codon_parts = []
    for codon, n in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])):
        aa = translate(codon, args.translation_table)
        reads_aa[aa] = reads_aa.get(aa, 0) + n
        codon_parts.append(f"{codon}:{n}:{aa}")
    reads_class = {c: 0 for c in CLASSES}
    for aa, n in reads_aa.items():
        reads_class[aa_class(aa, mutation)] += n
    depth = sum(reads_aa.values())

    flags = []
    locus_ploidy = 1 if "all" in args.haploid or mutation.chrom in args.haploid else 2
    norm_ploidy = 1 if "all" in args.haploid else 2
    genotype = "NA"
    if depth < args.min_depth:
        flags.append("low_depth")
    elif args.coverage and rcn is None:
        flags.append("rcn_NA")
    else:
        n_copies = total_copies(rcn, locus_ploidy, norm_ploidy, args.dup_rcn)
        copies = allocate_copies(reads_aa, n_copies, args.min_reads, args.min_frac)
        copies_class = {c: 0 for c in CLASSES}
        for aa, n in copies.items():
            copies_class[aa_class(aa, mutation)] += n
            if n == 0:
                flags.append(f"unplaced_{aa}")
        genotype = genotype_string(copies_class) or "NA"
        if reads_class["R"] and not copies_class["R"]:
            flags.append("R_below_threshold")
    if n_indel:
        flags.append(f"indel_reads_{n_indel}")

    rcn_str = f"{rcn:.4g}" if rcn is not None else "NA"
    return ";".join([
        genotype, str(reads_class["S"]), str(reads_class["R"]), str(reads_class["O"]),
        str(depth), rcn_str, ",".join(codon_parts) or "NA", ",".join(flags) or ".",
    ])


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run(args) -> None:
    """Entry point for ``ardu genotype``."""
    start_time = time.time()
    bam.configure_read_filters(args.min_mapq, args.min_baseq)

    for path in (args.bam, args.mutation, args.fasta, args.coverage):
        if path and not os.path.isfile(path):
            sys.exit(f"File not found: {path}")
    out_parent = os.path.dirname(args.outfile)
    if out_parent and not os.path.isdir(out_parent):
        sys.exit(f"Output directory does not exist: {out_parent}")

    with open(args.bam) as fh:
        bam_files = [line.strip() for line in fh if line.strip()]
    missing = [b for b in bam_files if not os.path.isfile(b)]
    if missing:
        sys.exit("BAM files not found:\n  " + "\n  ".join(missing))
    samples = [os.path.basename(b).removesuffix(".bam") for b in bam_files]

    try:
        mutations = read_mutations(args.mutation)
    except ValueError as exc:
        sys.exit(str(exc))
    if not mutations:
        sys.exit(f"No mutation in {args.mutation}.")

    # the reference codon must give the susceptible amino acid of the name,
    # otherwise a wrong coordinate or strand would give wrong genotypes silently
    ref_codons = {}
    errors = []
    with pysam.FastaFile(args.fasta) as fasta:
        for m in mutations:
            if m.chrom not in fasta.references:
                errors.append(f"{m.id}: chromosome {m.chrom} not in {args.fasta}")
                continue
            codon = reference_codon(fasta, m)
            aa = translate(codon, args.translation_table)
            if aa != m.ref_aa:
                errors.append(f"{m.id}: the reference codon {codon} gives {aa}, not {m.ref_aa}")
            ref_codons[m.id] = (codon, aa)
    if errors:
        sys.exit("Mutation file does not match the reference:\n  " + "\n  ".join(errors))

    rcn_table = {}
    if args.coverage:
        rcn_table = _read_rcn(args.coverage)
        absent = sorted({m.gene for m in mutations} - set(rcn_table))
        if absent:
            sys.exit(f"Genes not in {args.coverage}: {', '.join(absent)}")
        absent = sorted({s for m in mutations for s in samples if s not in rcn_table[m.gene]})
        if absent:
            sys.exit(f"Samples not in {args.coverage}: {', '.join(absent)}")

    cells = {m.id: {} for m in mutations}
    for bam_file, sample in zip(bam_files, samples):
        with pysam.AlignmentFile(bam_file, "rb") as samfile:
            for m in mutations:
                counts, n_indel = read_codons(samfile, m)
                rcn = rcn_table[m.gene][sample] if args.coverage else None
                cells[m.id][sample] = _genotype_cell(m, counts, n_indel, rcn, args)
        print(f"{sample}: {len(mutations)} codons read")

    out = f"{args.outfile}_codons.tsv"
    ploidy = ("all chromosomes haploid" if "all" in args.haploid else
              f"haploid: {', '.join(args.haploid)}; others and the normalisation diploid"
              if args.haploid else "all diploid")
    with open(out, "w") as fh:
        fh.write("# ArDu genotype outfile\n")
        fh.write(f"# Generated on: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")
        fh.write("# Format: Genotype;ReadsS;ReadsR;ReadsO;Depth;RCN;Codons;Flags\n")
        fh.write("# Genotype: copies of the susceptible (S, reference), resistant (R, named in the "
                 "mutation) and other (O) amino acids, e.g. 2S1R.\n")
        fh.write("# ReadsS, ReadsR, ReadsO: codons read for each class; Depth: their sum.\n")
        fh.write("# RCN:      relative copy number of the gene from the coverage file (NA without it).\n")
        fh.write("# Codons:   codon:reads:amino acid, on the coding strand.\n")
        fh.write("# Flags:    low_depth, rcn_NA, R_below_threshold, unplaced_<aa> (kept but no copy "
                 "left), indel_reads_<n> (reads with an indel in the codon, left out).\n")
        fh.write(f"# Read filters: MAPQ >= {args.min_mapq}, base quality >= {args.min_baseq}; "
                 f"translation table {args.translation_table}.\n")
        fh.write(f"# Copies: duplicated above {args.dup_rcn:g} copies per expected copy; an amino acid "
                 f"needs {args.min_reads} reads and max({args.min_frac:g}, half a copy) of them; "
                 f"genotype NA below {args.min_depth} codons. Ploidy: {ploidy}.\n")
        if not args.coverage:
            fh.write("# No coverage file: every locus taken as not duplicated.\n")
        for m in mutations:
            codon, aa = ref_codons[m.id]
            fh.write(f"# {m.id}: reference codon {codon} ({aa})\n")
        fh.write("mutation\t" + "\t".join(samples) + "\n")
        for m in mutations:
            fh.write(m.id + "\t" + "\t".join(cells[m.id][s] for s in samples) + "\n")
    print(f"Codon genotypes saved to {out}")
    print(f"Elapsed time: {round(time.time() - start_time, 1)} seconds")
