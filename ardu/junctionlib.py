from __future__ import annotations
import os
import ahocorasick
import pysam


def build_softclip_automaton(softclips_by_bkp, allowed_breakpoints=None):
    automaton = ahocorasick.Automaton()
    for bkp, clips in softclips_by_bkp.items():
        if allowed_breakpoints is not None and bkp not in allowed_breakpoints:
            continue
        for clip in clips:
            automaton.add_word(clip, (clip, bkp))
    automaton.make_automaton()
    return automaton


def find_reads_with_multiple_softclips_aho(bam_file, automaton, region=None, allowed_pairs=None):
    results = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        mapped_reads = list(bam.fetch(region[0], region[1], region[2])) if region else list(bam.fetch(until_eof=True))
        unmapped_reads = [r for r in bam.fetch(until_eof=True) if r.is_unmapped]
    for read in mapped_reads + unmapped_reads:
        if read.query_sequence is None:
            continue
        matches = []
        for end_pos, (pattern, bkp) in automaton.iter(read.query_sequence):
            start_pos = end_pos - len(pattern) + 1
            matches.append((start_pos, end_pos, bkp, pattern))
        for i in range(len(matches)):
            for j in range(i + 1, len(matches)):
                s1, e1, bkp1, seq1 = matches[i]
                s2, e2, bkp2, seq2 = matches[j]
                if bkp1 == bkp2:
                    continue
                if allowed_pairs is not None and (bkp1, bkp2) not in allowed_pairs:
                    continue
                seq_list = list(read.query_sequence.lower())
                for s, e, _, _ in [matches[i], matches[j]]:
                    for k in range(s, e + 1):
                        seq_list[k] = seq_list[k].upper()
                results.append({
                    "read_name": read.query_name,
                    "breakpoints": (bkp1, bkp2),
                    "capitalised_seq": "".join(seq_list),
                    "positions": [(s1, e1, bkp1, seq1), (s2, e2, bkp2, seq2)],
                })
                break
    return results


def write_junction_reads_fasta(junction_reads, output_dir, bam_name):
    pair_to_fh = {}
    for r in junction_reads:
        b1, b2 = r["breakpoints"]
        key = f"{b1}-{b2}"
        if key not in pair_to_fh:
            pair_to_fh[key] = open(os.path.join(output_dir, f"{bam_name}_{key}_junction_reads.fasta"), "w")
        pair_to_fh[key].write(f">{bam_name}:{r['read_name']}:{b1}-{b2}\n{r['capitalised_seq']}\n")
    for fh in pair_to_fh.values():
        fh.close()


def write_junction_reads_fastq(bam_file, junction_reads, output_dir, bam_name):
    pair_to_fh = {}
    read_name_set = {r["read_name"] for r in junction_reads}
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.query_name not in read_name_set:
                continue
            quals = read.query_qualities or [30] * len(read.query_sequence)
            qual_str = "".join(chr(q + 33) for q in quals)
            for r in junction_reads:
                if r["read_name"] != read.query_name:
                    continue
                b1, b2 = r["breakpoints"]
                key = f"{b1}-{b2}"
                if key not in pair_to_fh:
                    pair_to_fh[key] = open(os.path.join(output_dir, f"{bam_name}_{key}_junction_reads.fastq"), "w")
                pair_to_fh[key].write(f"@{bam_name}:{read.query_name}:{b1}-{b2}\n{read.query_sequence}\n+\n{qual_str}\n")
                break
    for fh in pair_to_fh.values():
        fh.close()
