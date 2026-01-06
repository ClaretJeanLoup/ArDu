import pysam
import argparse
import os
import re
from collections import defaultdict
from Bio.Seq import Seq  # translate amino acids

def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

def translate_codon(codon, table=1):
    try:
        return str(Seq(codon.upper()).translate(table=table))
    except Exception as e:
        print(f"[WARNING] Could not translate codon {codon}: {e}")
        return "X"

def extract_codon_from_pileup(bam_file, chrom, start, stop, strand, table=1, min_mapq=40, min_base_quality=20):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    codons = defaultdict(int)
    positions = [start - 1, start, start + 1]
    read_bases = defaultdict(dict)

    for pileupcolumn in samfile.pileup(chrom, positions[0], positions[-1] + 1, truncate=True):
        ref_pos = pileupcolumn.reference_pos
        if ref_pos not in positions:
            continue

        for pileupread in pileupcolumn.pileups:
            if pileupread.is_del or pileupread.is_refskip:
                continue
            read = pileupread.alignment
            read_pos = pileupread.query_position
            if read.mapping_quality < min_mapq or read_pos is None:
                continue
            base_qual = read.query_qualities[read_pos]
            if base_qual < min_base_quality:
                continue
            base = read.query_sequence[read_pos]
            read_bases[read.query_name][ref_pos] = base

    for read_name, base_dict in read_bases.items():
        if all(pos in base_dict for pos in positions):
            codon = ''.join(base_dict[pos] for pos in positions)
            if strand == "reverse":
                codon = reverse_complement(codon)
            codons[codon] += 1

    samfile.close()
    return codons

def extract_reference_codon(reference_fasta, chrom, start, stop, strand):
    ref = pysam.FastaFile(reference_fasta)
    ref_seq = ref.fetch(chrom, start - 1, stop)
    ref.close()
    return reverse_complement(ref_seq) if strand == "reverse" else ref_seq

def parse_coverage_file(coverage_file, bam_files):
    bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in bam_files]

    with open(coverage_file) as f:
        header = None
        while header is None or header[0].startswith('#'):
            header = f.readline().strip().split('\t')
        col_indices = {col.strip(): idx for idx, col in enumerate(header)}

        coverage_data = defaultdict(dict)

        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            gene = parts[0].strip()
            for sample in bam_basenames:
                if sample not in col_indices:
                    continue

                idx = col_indices[sample]
                values = parts[idx].strip().split(';')

                try:
                    copy_number = float(values[-1])  # Take the last value
                    coverage_data[gene][sample] = copy_number
                except Exception as e:
                    print(f"[WARNING] Could not parse copy number for {gene}/{sample}: {values[-1]}")
                    coverage_data[gene][sample] = None

    return coverage_data

def parse_mutation_name(mutation_name):
    # Supports both single mutations and multiple mutations separated by /
    match = re.match(r"([A-Z])(\d+)([A-Z/]+)", mutation_name)
    if not match:
        raise ValueError(f"Invalid mutation format: {mutation_name}")
    ref_aa = match.group(1)
    mutation_aas = set(match.group(3).split('/')) 
    return ref_aa, mutation_aas

def determine_pcp(ref_aa, mutation_aas, observed_aas, copy_number, threshold=1.4):
    # If observed amino acids have neither S nor R
    if ref_aa not in observed_aas and not (observed_aas & mutation_aas):
        return "NA"  

    if ref_aa in observed_aas and not (observed_aas & mutation_aas):
        if copy_number > threshold:
            return "Sx"
        return "SS"  # Only S

    # If R is observed and ref_aa is not observed
    if not ref_aa in observed_aas and (observed_aas & mutation_aas):
        if copy_number > threshold:
            return "Rx"  # Duplicated resistance
        else:
            return "RR"  # Resistance only

    # If both ref_aa and mutation_aas are observed D or mixed R and S alleels
    if ref_aa in observed_aas and (observed_aas & mutation_aas):
        if copy_number > threshold:
            return "D*"  # Can't predict so set as D*
        else:
            return "RS" 

    return "NA"

def parse_codon_string(codon_string):
    """
    Parses codon output string into codons and fields
    """
    codons = []
    fields = {}
    for part in codon_string.split(';'):
        if part.startswith("TOT:") or part.startswith("REF:") or part.startswith("Predicted:") or part.startswith("CopyNumber:") or part.startswith("PCP:") or part.startswith("TotS:") or part.startswith("TotR:"):
            key, value = part.split(':', 1)
            fields[key] = value
        elif part: 
            codons.append(part)
    return codons, fields

def main():
    parser = argparse.ArgumentParser(description="Extract codons from BAM files.")
    parser.add_argument("--bam", required=True, help="Text file listing BAM files, one per line")
    parser.add_argument("--mutation", required=True, help="Tab-separated file with chromosome, start, stop, gene, mutation name, strand (forward/reverse)")
    parser.add_argument("--outfile", required=True, help="Output file prefix")
    parser.add_argument("--reference", help="Reference genome in FASTA format")
    parser.add_argument("--coverage", help="Coverage file for normalisation")
    parser.add_argument("--translation-table", type=int, default=1, help="NCBI translation table number (default: 1 = Standard Code)")
    parser.add_argument("--pcp", action="store_true", help="Whether to calculate Partialy Constructed Phenotype (PCP)")
    parser.add_argument("--tab", action="store_true", help="Output in tabular format with columns per BAM")
    parser.add_argument("--min-mapq", type=int, default=40, help="Minimum mapping quality (default: 40)")
    parser.add_argument("--min-base-quality", type=int, default=20, help="Minimum base quality (default: 20)")

    args = parser.parse_args()

    with open(args.bam) as bf:
        bam_files = bf.read().splitlines()
        bam_basenames = [os.path.splitext(os.path.basename(bam))[0] for bam in bam_files]

    coverage_data = {}
    if args.coverage:
        coverage_data = parse_coverage_file(args.coverage, bam_files)

    mutation_data = {}
    with open(args.mutation) as mf:
        for line in mf:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                print(f"[WARNING] Incomplete line skipped: {line}")
                continue

            chrom, start, stop, gene, mutation_name, strand = parts[:6]
            strand = strand.strip().lower()
            if strand not in ("forward", "reverse"):
                print(f"[WARNING] Invalid strand '{strand}' in line: {line}")
                strand = "forward"

            start, stop = int(start), int(stop)
            mutation_id = f"{mutation_name}:{gene}:{chrom}:{start}-{stop}:{strand}"

            if mutation_id not in mutation_data:
                mutation_data[mutation_id] = {}

            ref_aa, mutation_aas = parse_mutation_name(mutation_name)

            for bam_file, bam_name in zip(bam_files, bam_basenames):
                try:
                    codons = extract_codon_from_pileup(
                        bam_file, chrom, start, stop, strand,
                        table=args.translation_table,
                        min_mapq=args.min_mapq,
                        min_base_quality=args.min_base_quality)
                    total = sum(codons.values())

                    codon_entries = []
                    aa_counts = defaultdict(int)
                    observed_aas = set()
                    for codon, count in codons.items():
                        aa = translate_codon(codon, table=args.translation_table)
                        aa_counts[aa] += count
                        observed_aas.add(aa)
                        codon_entries.append(f"{codon}:{count}:{aa}")

                    codon_str = ";".join(codon_entries) + f";TOT:{total}"

                    if args.reference:
                        ref_codon = extract_reference_codon(args.reference, chrom, start, stop, strand)
                        ref_aa_seq = translate_codon(ref_codon, table=args.translation_table)
                        codon_str += f";REF:{ref_codon}:{ref_aa_seq}"

                    if gene in coverage_data and bam_name in coverage_data[gene] and coverage_data[gene][bam_name] is not None:
                        copy_number = coverage_data[gene][bam_name]
                        if total > 0:
                            frac_counts = {
                                aa: (count / total) * (2 * copy_number)
                                for aa, count in aa_counts.items()}

                            predicted = "/".join([f"{round(frac_counts[aa])}{aa}" for aa in sorted(frac_counts.keys())])
                            codon_str += f";Predicted:{predicted};CopyNumber:{copy_number:.2f}"

                            # PCP uses known R codons
                            relevant_aas = {aa for aa in frac_counts if aa == ref_aa or aa in mutation_aas}
                            if args.pcp:
                                pcp = determine_pcp(ref_aa, mutation_aas, relevant_aas, copy_number)
                                codon_str += f";PCP:{pcp}"

                                tot_s = frac_counts.get(ref_aa, 0)
                                tot_r = sum(frac_counts[aa] for aa in mutation_aas if aa in frac_counts)
                                tot_alt = sum(frac_counts[aa] for aa in frac_counts if aa != ref_aa and aa not in mutation_aas) #alt that are not R 
                                codon_str += f";TotAlt:{tot_alt:.2f};TotS:{tot_s:.2f};TotR:{tot_r:.2f}"
                    else:
                        codon_str += ";Predicted:NA;CopyNumber:NA;PCP:NA;TotAlt:NA;TotS:NA;TotR:NA"


                    mutation_data[mutation_id][bam_name] = codon_str

                except Exception as e:
                    print(f"[ERROR] Processing {bam_file} failed: {e}")
                    mutation_data[mutation_id][bam_name] = "ERROR"

    out_file = f"{args.outfile}_codons.tsv"
    with open(out_file, 'w') as outf:
        if args.tab:
            headers = ["mutation_id"]
            for bam in bam_basenames:
                headers += [f"codons_{bam}", f"TOT_{bam}", f"REF_{bam}", f"Predicted_{bam}", f"CopyNumber_{bam}", f"PCP_{bam}", f"TotS_{bam}", f"TotR_{bam}"]
            outf.write("\t".join(headers) + "\n")

            for mutation_id, bam_vals in mutation_data.items():
                row = [mutation_id]
                for bam in bam_basenames:
                    codon_string = bam_vals.get(bam, "NA")
                    if codon_string == "NA" or codon_string == "ERROR":
                        row += ["NA"] * 8  # 8 fields : codons + 7 
                        continue
                    codons, fields = parse_codon_string(codon_string)
                    codon_joined = ";".join(codons) if codons else "NA"
                    row.append(codon_joined)
                    for key in ["TOT", "REF", "Predicted", "CopyNumber", "PCP", "TotS", "TotR"]:
                        row.append(fields.get(key, "NA"))
                outf.write("\t".join(row) + "\n")
        else:
            # Regular output mode
            outf.write("mutation_id\t" + "\t".join(bam_basenames) + "\n")
            for mutation_id, bam_vals in mutation_data.items():
                row = [mutation_id] + [bam_vals.get(bam, "NA") for bam in bam_basenames]
                outf.write("\t".join(row) + "\n")

    print(f"Codon data saved to {out_file}")


if __name__ == "__main__":
    main()
