import pysam
import argparse
import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import ahocorasick

def extract_breakpoint_from_qseqid(qseqid):
    match = re.search(r'_(\d+)_', qseqid)
    return int(match.group(1)) if match else None

def filter_blast_hits_by_breakpoint(blast_output, breakpoints_df, extension, filtered_output, pairs_file=None):
    allowed_pairs = set()
    if pairs_file:
        pairs_df = pd.read_csv(pairs_file, sep="\t")
        for _, row in pairs_df.iterrows():
            source = int(row["source_bkp"])
            target = int(row["target_bkp"])
            allowed_pairs.add((source, target))
            allowed_pairs.add((target, source))

    regions_by_bkp = {}
    for _, row in breakpoints_df.iterrows():
        chrom = row["chromosome"]
        pos = row["position"]
        bkp = row["breakpoint_number"]
        regions_by_bkp[bkp] = (chrom, pos - extension, pos + extension)

    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
    blast_df = pd.read_csv(blast_output, sep="\t", names=columns)

    filtered_hits = []
    for _, row in blast_df.iterrows():
        qseqid = row["qseqid"]
        sseqid = str(row["sseqid"])
        sstart = min(row["sstart"], row["send"])
        send = max(row["sstart"], row["send"])

        source_bkp = extract_breakpoint_from_qseqid(qseqid)
        if source_bkp is None:
            continue

        for target_bkp, (chrom, start, end) in regions_by_bkp.items():
            if target_bkp == source_bkp:
                continue
            if pairs_file and (source_bkp, target_bkp) not in allowed_pairs:
                continue
            if sseqid == chrom and not (send < start or sstart > end):
                filtered_hits.append(row)
                break

    pd.DataFrame(filtered_hits).to_csv(filtered_output, sep="\t", index=False, header=False)
    print(f"Filtered BLAST results (with allowed pairs) saved to {filtered_output}")

def extract_soft_clipped_positions(bam_file, region, min_clip, output_fasta, output_tsv, breakpoint):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        with open(output_fasta, "a") as fasta_out, open(output_tsv, "a") as tsv_out:
            for read in bam.fetch(region[0], region[1], region[2]):
                if read.cigar is None:
                    continue

                mapping_position = read.reference_start + 1
                cigar = read.cigarstring
                chrom = read.reference_name

                mate_position = read.next_reference_start if not read.mate_is_unmapped else None
                insert_size = abs(mapping_position - mate_position) if mate_position is not None else None

                if read.cigar[0][0] == 4 and read.cigar[0][1] >= min_clip:
                    left_clip = read.query_sequence[:read.cigar[0][1]]
                    fasta_out.write(f">{chrom}:{mapping_position}_{read.query_name}_{breakpoint}_L\n{left_clip}\n")
                    tsv_out.write(f"{bam_file}\t{read.query_name}\t{mapping_position}\t{mapping_position}\tleftmost_clip\t{cigar}\t{mate_position}\t{insert_size}\t{region}\t{breakpoint}\n")

                if read.cigar[-1][0] == 4 and read.cigar[-1][1] >= min_clip:
                    right_clip = read.query_sequence[-read.cigar[-1][1]:]
                    clip_position = mapping_position + (read.query_length - read.cigar[-1][1])
                    fasta_out.write(f">{chrom}:{clip_position}_{read.query_name}_{breakpoint}_R\n{right_clip}\n")
                    tsv_out.write(f"{bam_file}\t{read.query_name}\t{mapping_position}\t{clip_position}\trightmost_clip\t{cigar}\t{mate_position}\t{insert_size}\t{region}\t{breakpoint}\n")
    bam.close()

def run_blast(query_fasta, reference, blast_output):
    if not os.path.exists(f"{reference}.nhr"):
        print("Creating BLAST database...")
        subprocess.run(["makeblastdb", "-in", reference, "-dbtype", "nucl"], check=True)

    print("Running BLAST...")
    subprocess.run(["blastn", "-query", query_fasta, "-db", reference, "-out", blast_output, "-outfmt", "6"], check=True)
    print(f"BLAST results saved in {blast_output}")

def build_softclip_automaton(softclips_by_bkp):
    A = ahocorasick.Automaton()
    for bkp, clips in softclips_by_bkp.items():
        for clip in clips:
            A.add_word(clip, (clip, bkp))  # store pattern + bkp
    A.make_automaton()
    return A

def find_reads_with_multiple_softclips_aho(bam_file, automaton, region=None):
    bam = pysam.AlignmentFile(bam_file, "rb")
    results = []

    if region:
        chrom, start, end = region
        mapped_reads = bam.fetch(chrom, start, end)
    else:
        mapped_reads = bam.fetch(until_eof=True)

    unmapped_reads = [read for read in bam.fetch(until_eof=True) if read.is_unmapped]
    for read in list(mapped_reads) + unmapped_reads:
        if read.query_sequence is None:
            continue

        matches = []
        for end_pos, (pattern, bkp) in automaton.iter(read.query_sequence):
            start_pos = end_pos - len(pattern) + 1
            matches.append((start_pos, end_pos, bkp))

        # fnd any valid pair (withoverlap)
        for i in range(len(matches)):
            for j in range(i + 1, len(matches)):
                s1, e1, bkp1 = matches[i]
                s2, e2, bkp2 = matches[j]
                if bkp1 == bkp2:
                    continue
                capitalised = capitalise_matches(read.query_sequence, [(s1, e1, bkp1), (s2, e2, bkp2)])
                # get pair in read order
                if s1 < s2:
                    bkps_in_order = [bkp1, bkp2]
                else:
                    bkps_in_order = [bkp2, bkp1]
                results.append((read.query_name, bkps_in_order, capitalised))
                break  # one valid pair is enough

    bam.close()
    return results


def capitalise_matches(seq, matches):
    """
    seq: original read sequence
    matches: list of (start, end, bkp)
    Returns: string where matched intervals are uppercase, others are lowercase
    """
    seq = seq.lower()
    seq_list = list(seq)
    for s, e, _ in matches:
        for i in range(s, e + 1):
            seq_list[i] = seq_list[i].upper()
    return "".join(seq_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract soft-clipped sequences and scan for junction reads.")
    parser.add_argument("-b", "--bam_list", required=True, help="File listing BAM files.")
    parser.add_argument("-i", "--input", required=True, help="Breakpoint file.")
    parser.add_argument("-o", "--output", required=True, help="Output prefix.")
    parser.add_argument("-s", "--size", type=int, default=30, help="Minimum soft-clip length.")
    parser.add_argument("-e", "--extension", type=int, default=30, help="Region around breakpoint.")
    parser.add_argument("--blast", metavar="REFERENCE", help="Path to reference genome for BLAST.")
    parser.add_argument("--pairs", help="TSV with allowed breakpoint pairs.")
    parser.add_argument("--junction", action="store_true", help="Search for junction reads.")
    parser.add_argument("--junction-region", help="Genomic region (chr:start-end) for junction search.")

    args = parser.parse_args()
    fasta_output = f"{args.output}_soft_clipped_sequences.fasta"
    blast_output = f"{args.output}_BLAST_results.txt"
    tsv_output = f'{args.output}.tsv'

    bam_mapping = {}
    with open(args.bam_list) as f:
        for line in f:
            full_path = line.strip()
            base = os.path.basename(full_path)
            if base.endswith(".bam"):
                bam_mapping[base] = full_path

    breakpoints_df = pd.read_csv(args.input, sep="\t")
    breakpoints_df["bam_file"] = breakpoints_df["bam_file"].apply(lambda x: os.path.basename(x) + ".bam")

    with open(tsv_output, "w") as out:
        out.write("bam_file\tRead_name\tPosition\tClipping_position\tClipping_end\tCIGAR\tMate_position\tInsert_size\tRegion\tBreakpoint\n")

    softclips_by_bkp = {}

    for _, row in breakpoints_df.iterrows():
        bam_basename = row["bam_file"]
        if bam_basename not in bam_mapping:
            raise ValueError(f"BAM file {bam_basename} not found in --bam_list.")

        bam_file = bam_mapping[bam_basename]
        chrom = row["chromosome"]
        pos = row["position"]
        bkp = row["breakpoint_number"]
        region = (chrom, pos - args.extension, pos + args.extension)

        extract_soft_clipped_positions(bam_file, region, args.size, fasta_output, tsv_output, bkp)

        if bkp not in softclips_by_bkp:
            softclips_by_bkp[bkp] = []
        bam = pysam.AlignmentFile(bam_file, "rb")
        for read in bam.fetch(region[0], region[1], region[2]):
            if read.cigar is None or read.query_sequence is None:
                continue
            if read.cigar[0][0] == 4 and read.cigar[0][1] >= args.size:
                softclips_by_bkp[bkp].append(read.query_sequence[:read.cigar[0][1]])
            if read.cigar[-1][0] == 4 and read.cigar[-1][1] >= args.size:
                softclips_by_bkp[bkp].append(read.query_sequence[-read.cigar[-1][1]:])
        bam.close()

    if args.blast:
        run_blast(fasta_output, args.blast, blast_output)
        filtered_output = f"{args.output}_BLAST_filtered_pairs.txt"
        filter_blast_hits_by_breakpoint(blast_output, breakpoints_df, args.extension, filtered_output, args.pairs)

    if args.junction:
        print("Running Aho-Corasick junction scan...")
        automaton = build_softclip_automaton(softclips_by_bkp)

        junction_region = None
        if args.junction_region:
            try:
                chrom, start_end = args.junction_region.split(":")
                start, end = map(int, start_end.split("-"))
                junction_region = (chrom, start, end)
            except ValueError:
                raise ValueError("Invalid format for --junction-region. Use chr:start-end (e.g., chr1:100000-200000)")

        with open(args.bam_list) as f:
            bam_files = [line.strip() for line in f if line.strip()]

        junction_output = f"{args.output}_junction_reads.fasta"
        with open(junction_output, "w") as out:
            for bam in bam_files:
                matches = find_reads_with_multiple_softclips_aho(bam, automaton, junction_region)
                for read_name, bkps, capitalised_seq in matches:
                    out.write(f">{os.path.basename(bam)}:{read_name}:{'-'.join(map(str, bkps))}\n{capitalised_seq}\n")


        print(f"Junction reads saved to {junction_output}")
