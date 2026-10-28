"""
``ardu regions`` — build the ``-r`` region file of ``ardu coverage`` from gene names
and an annotation, or from query sequences blasted against the reference.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pysam

from ardu.locate import (
    Annotation,
    alternative_hits,
    best_two_chains,
    hsps_from_blast,
    merge_intervals,
    parse_region_line,
)

BLAST_FIELDS = "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen btop"


# ---------------------------------------------------------------------------
# Argparse registration
# ---------------------------------------------------------------------------

def register(subparsers) -> None:
    """Add the ``regions`` sub-command to *subparsers*."""
    parser = subparsers.add_parser(
        "regions",
        help="Build the -r region file from gene names and a gff/gtf, or by blast.",
        description=(
            "ArDu regions — write the 4-column region file of ardu coverage "
            "(chromosome, start, stop, locus_name; 1-based, inclusive), either from a "
            "list of gene names and a gff3/gtf annotation, or from query sequences "
            "blasted against the reference."
        ),
    )

    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("-g", "--gff", metavar="ANNOTATION",
                        help="gff3 or gtf annotation (may be gzipped). Needs -l.")
    source.add_argument("-q", "--query", metavar="FASTA",
                        help="Multifasta of locus sequences (cDNA, CDS or genomic), one per "
                             "locus, named by the first word of the header. Needs --reference.")
    parser.add_argument("-l", "--loci", metavar="LIST",
                        help="With -g: one gene or transcript per line, optionally followed by "
                             "the locus name to use (lines sharing a locus name are pooled). "
                             "A line may also be a region, 'chrom:start-end locus_name'.")
    parser.add_argument("--reference", metavar="FASTA",
                        help="With -q: reference fasta of the bams. A blast database is built "
                             "next to it if missing.")
    parser.add_argument("-o", "--outfile", required=True, metavar="PREFIX",
                        help="Prefix of the output files: <prefix>_regions.bed, and with -q "
                             "<prefix>_blast.tsv and <prefix>_hits.tsv.")

    parser.add_argument("--feature", choices=["exon", "CDS", "gene"], default="exon",
                        help="Part of each locus to keep: its exons or CDS (merged across "
                             "isoforms), or the whole gene span. With -q, 'exon' and 'CDS' "
                             "keep the blast hits and 'gene' their whole span.")
    parser.add_argument("--pad", type=int, default=0,
                        help="Bases added on each side of every interval.")
    parser.add_argument("--chrom-map", metavar="FILE",
                        help="Two columns: chromosome name in the annotation, name in the bams "
                             "(e.g. NT_078266.2 AgamP4_2R). Unlisted names are kept.")
    parser.add_argument("--fasta", metavar="REFERENCE",
                        help="Reference fasta (indexed) of the bams: checks the chromosome names "
                             "and clips padded intervals to the chromosome ends. With -q the "
                             "--reference is used.")
    parser.add_argument("--plot-margin", type=int, metavar="BP",
                        help="Also write <prefix>_plot.bed for --plot-interval, each locus span "
                             "extended by this many bases on each side.")

    blast = parser.add_argument_group("blast (-q)")
    blast.add_argument("--task", default="megablast",
                       choices=["megablast", "dc-megablast", "blastn"],
                       help="blastn task; dc-megablast or blastn for sequences of another species.")
    blast.add_argument("--evalue", type=float, default=1e-10)
    blast.add_argument("--threads", type=int, default=1)
    blast.add_argument("--max-intron", type=int, default=50000,
                       help="Largest gap (bp) on the reference between two hits of one chain.")
    blast.add_argument("--min-intron", type=int, default=20,
                       help="Hits are cut where the reference has a run of at least this many "
                            "bases missing from the query (an intron bridged by a gapped hit).")
    blast.add_argument("--hit-overlap", type=int, default=30,
                       help="Overlap (bp) allowed between consecutive hits, on the query and "
                            "on the reference.")
    blast.add_argument("--min-cover", type=float, default=0.8,
                       help="Warn when the best chain covers less of the query than this.")
    blast.add_argument("--second-ratio", type=float, default=0.9,
                       help="Warn when a second chain elsewhere, or a hit elsewhere for part of "
                            "the chain, scores at least this fraction of the chosen one "
                            "(another copy in the reference).")

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def _read_list(path: str) -> list[tuple[str, str]]:
    """``(entry, locus_name)`` pairs of the -l file."""
    entries = []
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if not fields or fields[0].startswith("#"):
                continue
            entries.append((fields[0], fields[1] if len(fields) > 1 else fields[0]))
    return entries


def _read_chrom_map(path: str | None) -> dict[str, str]:
    if not path:
        return {}
    out = {}
    with open(path) as fh:
        for line in fh:
            fields = line.split()
            if len(fields) >= 2 and not fields[0].startswith("#"):
                out[fields[0]] = fields[1]
    return out


def _query_names(path: str) -> list[str]:
    names = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        sys.exit(f"Query names must be unique: {', '.join(dup)}")
    return names


# ---------------------------------------------------------------------------
# Annotation mode
# ---------------------------------------------------------------------------

def _from_annotation(args) -> dict[str, list[tuple[str, int, int]]]:
    if not args.loci:
        sys.exit("-g needs -l, the list of genes.")
    entries = _read_list(args.loci)
    if not entries:
        sys.exit(f"No gene in {args.loci}.")
    genes = [e for e, _ in entries if parse_region_line(e) is None]
    ann = Annotation(args.gff) if genes else None

    loci: dict[str, list[tuple[str, int, int]]] = {}
    errors = []
    for entry, locus in entries:
        region = parse_region_line(entry)
        if region is not None:
            loci.setdefault(locus, []).append(region)
            continue
        found = ann.find(entry)
        if not found:
            errors.append(f"{entry}: not found in {args.gff}.")
            continue
        if len(found) > 1:
            errors.append(f"{entry}: matches {len(found)} features, give one of their IDs: "
                          + "; ".join(f.label() for f in found))
            continue
        intervals = ann.intervals(found[0], args.feature)
        if not intervals:
            errors.append(f"{entry}: no {args.feature} under {found[0].label()}; "
                          "try --feature gene.")
            continue
        print(f"{entry} -> {found[0].label()}, {len(intervals)} interval(s) as {locus}")
        loci.setdefault(locus, []).extend(intervals)
    if errors:
        sys.exit("\n".join(errors))
    return loci


# ---------------------------------------------------------------------------
# Blast mode
# ---------------------------------------------------------------------------

def _run_blast(args, blast_output: str) -> None:
    ref = args.reference
    if not any(os.path.exists(f"{ref}.{ext}") for ext in ("nhr", "nal", "ndb")):
        print("Creating BLAST database...")
        subprocess.run(["makeblastdb", "-in", ref, "-dbtype", "nucl"], check=True)
    print("Running BLAST...")
    subprocess.run(
        ["blastn", "-task", args.task, "-query", args.query, "-db", ref,
         "-evalue", str(args.evalue), "-num_threads", str(args.threads),
         "-outfmt", f"6 {BLAST_FIELDS}", "-out", blast_output],
        check=True,
    )


def _from_blast(args, prefix: str) -> dict[str, list[tuple[str, int, int]]]:
    if not args.reference:
        sys.exit("-q needs --reference, the fasta to blast against.")
    if args.loci:
        sys.exit("-l is for -g; with -q the loci are the query sequences.")
    names = _query_names(args.query)
    if not names:
        sys.exit(f"No sequence in {args.query}.")

    blast_output = f"{prefix}_blast.tsv"
    _run_blast(args, blast_output)

    hsps: dict[str, list] = {n: [] for n in names}
    qlen: dict[str, int] = {}
    with open(blast_output) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            hsps[fields[0]].extend(hsps_from_blast(fields, args.min_intron))
            qlen[fields[0]] = int(fields[10])

    loci: dict[str, list[tuple[str, int, int]]] = {}
    rows = []
    for name in names:
        first, second = best_two_chains(hsps[name], args.max_intron, args.hit_overlap)
        if first is None:
            print(f"WARNING: {name} has no hit; left out of the region file.")
            rows.append([name] + ["NA"] * 12)
            continue
        cover = first.query_cover(qlen[name])
        ratio = second.score / first.score if second else 0.0
        alts = alternative_hits(first, hsps[name], args.second_ratio)
        flags = []
        if cover < args.min_cover:
            flags.append("low_cover")
            print(f"WARNING: {name}: the best chain covers {cover:.0%} of the query.")
        if ratio >= args.second_ratio or alts:
            flags.append("second_copy")
        if ratio >= args.second_ratio:
            print(f"WARNING: {name}: a second chain at {second.chrom}:{second.start}-"
                  f"{second.end} scores {ratio:.0%} of the best one; check which copy you want.")
        for h, g in alts:
            print(f"WARNING: {name}: the hit at {h.chrom}:{h.sstart}-{h.send} also hits "
                  f"{g.chrom}:{g.sstart}-{g.send} ({g.bits / h.bits:.0%} of its bitscore): "
                  "another copy of this part of the locus.")
        if args.feature == "gene":
            intervals = [(first.chrom, first.start, first.end)]
        else:
            intervals = merge_intervals([(h.chrom, h.sstart, h.send) for h in first.hsps])
        loci[name] = intervals
        print(f"{name} -> {first.chrom}:{first.start}-{first.end} ({first.strand}), "
              f"{len(first.hsps)} hit(s), {cover:.0%} of the query")
        rows.append([
            name, first.chrom, first.strand, first.start, first.end, len(first.hsps),
            f"{cover:.3f}", f"{first.pident():.2f}", f"{first.score:g}",
            f"{second.chrom}:{second.start}-{second.end}" if second else "NA",
            f"{ratio:.3f}", len(alts), ",".join(flags) or ".",
        ])

    hits_output = f"{prefix}_hits.tsv"
    with open(hits_output, "w") as fh:
        fh.write("locus\tchrom\tstrand\tstart\tend\tn_hits\tquery_cover\tpident\tbitscore\t"
                 "second_chain\tsecond_ratio\talt_hits\tflags\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
    print(f"Chosen hits saved to {hits_output}")
    return loci


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def _finish(loci, chrom_map, pad, fasta_path):
    """Rename chromosomes, pad, clip to the chromosome ends and merge each locus."""
    lengths = None
    if fasta_path:
        with pysam.FastaFile(fasta_path) as fasta:
            lengths = dict(zip(fasta.references, fasta.lengths))
    out, missing = {}, set()
    for locus, intervals in loci.items():
        done = []
        for chrom, start, end in intervals:
            chrom = chrom_map.get(chrom, chrom)
            start, end = max(1, start - pad), end + pad
            if lengths is not None:
                if chrom not in lengths:
                    missing.add(chrom)
                    continue
                end = min(end, lengths[chrom])
            done.append((chrom, start, end))
        out[locus] = merge_intervals(done)
    if missing:
        sys.exit(f"Chromosomes not in {fasta_path}: {', '.join(sorted(missing))}. "
                 "Map them to the fasta names with --chrom-map.")
    return out


def run(args) -> None:
    prefix = args.outfile
    out_dir = os.path.dirname(os.path.abspath(prefix))
    if not os.path.isdir(out_dir):
        sys.exit(f"Output directory does not exist: {out_dir}")
    for path in (args.gff, args.query, args.loci, args.reference, args.fasta, args.chrom_map):
        if path and not os.path.exists(path):
            sys.exit(f"File not found: {path}")

    if args.gff:
        loci = _from_annotation(args)
        fasta = args.fasta
    else:
        loci = _from_blast(args, prefix)
        fasta = args.reference
    loci = _finish(loci, _read_chrom_map(args.chrom_map), args.pad, fasta)
    if not loci:
        sys.exit("No locus to write.")

    regions_output = f"{prefix}_regions.bed"
    with open(regions_output, "w") as fh:
        for locus, intervals in loci.items():
            for chrom, start, end in intervals:
                fh.write(f"{chrom}\t{start}\t{end}\t{locus}\n")
    print(f"Region file ({len(loci)} loci) saved to {regions_output}")

    if args.plot_margin is not None:
        plot_output = f"{prefix}_plot.bed"
        with open(plot_output, "w") as fh:
            for locus, intervals in loci.items():
                chroms = {c for c, _, _ in intervals}
                if len(chroms) > 1:
                    print(f"WARNING: {locus} spans several chromosomes; no plot interval.")
                    continue
                start = max(1, min(s for _, s, _ in intervals) - args.plot_margin)
                end = max(e for _, _, e in intervals) + args.plot_margin
                fh.write(f"{locus}\t{chroms.pop()}:{start}-{end}\n")
        print(f"Plot intervals saved to {plot_output}")
