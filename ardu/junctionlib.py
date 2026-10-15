from __future__ import annotations
import os
import ahocorasick
import pysam

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def build_softclip_automaton(softclips_by_bkp, allowed_breakpoints=None):
    automaton = ahocorasick.Automaton()
    for bkp, clips in softclips_by_bkp.items():
        if allowed_breakpoints is not None and bkp not in allowed_breakpoints:
            continue
        for clip in clips:
            automaton.add_word(clip, (clip, bkp))
    automaton.make_automaton()
    return automaton


def _revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def _junction_pairs(seq, automaton, allowed_pairs):
    """Every distinct pair of softclip matches from different breakpoints in *seq*."""
    matches = []
    for end_pos, (pattern, bkp) in automaton.iter(seq):
        matches.append((end_pos - len(pattern) + 1, end_pos, bkp, pattern))
    pairs = {}
    for i in range(len(matches)):
        for j in range(i + 1, len(matches)):
            bkp1, bkp2 = matches[i][2], matches[j][2]
            if bkp1 == bkp2 or (bkp1, bkp2) in pairs:
                continue
            if allowed_pairs is not None and (bkp1, bkp2) not in allowed_pairs:
                continue
            pairs[(bkp1, bkp2)] = (matches[i], matches[j])
    return pairs


def _scan_read(read, automaton, allowed_pairs, results):
    seq = read.query_sequence
    if seq is None:
        return
    quals = read.query_qualities
    # Softclip patterns come from mapped reads, stored in reference orientation.
    # Unmapped reads keep their sequencing orientation, so both strands are searched.
    candidates = [(seq, quals, "+")]
    if read.is_unmapped:
        candidates.append((_revcomp(seq), quals[::-1] if quals is not None else None, "-"))
    for s, q, strand in candidates:
        pairs = _junction_pairs(s, automaton, allowed_pairs)
        for (bkp1, bkp2), (m1, m2) in pairs.items():
            seq_list = list(s.lower())
            for start, end, _, _ in (m1, m2):
                for k in range(start, end + 1):
                    seq_list[k] = seq_list[k].upper()
            results.append({
                "read_name": read.query_name,
                "breakpoints": (bkp1, bkp2),
                "capitalised_seq": "".join(seq_list),
                "positions": [m1, m2],
                "strand": strand,
                "sequence": s,
                "qualities": q,
            })
        if pairs:
            return


def find_reads_with_multiple_softclips_aho(bam_file, automaton, region=None, allowed_pairs=None,
                                           include_unplaced=True):
    """Reads containing softclip sequences from two different breakpoints.

    With *region*, reads are streamed from that region only, plus (when
    *include_unplaced*) the unplaced unmapped reads, fetched through the index
    with ``"*"``. Without *region*, the whole BAM is streamed once.
    """
    results = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        if region:
            for read in bam.fetch(region[0], region[1], region[2]):
                _scan_read(read, automaton, allowed_pairs, results)
            if include_unplaced:
                for read in bam.fetch("*"):
                    _scan_read(read, automaton, allowed_pairs, results)
        else:
            for read in bam.fetch(until_eof=True):
                _scan_read(read, automaton, allowed_pairs, results)
    return results


def pair_label(b1, b2):
    """``1-2`` for plain breakpoint numbers, ``locus_1-2`` for ``(locus, number)`` keys."""
    if isinstance(b1, tuple):
        return f"{b1[0]}_{b1[1]}-{b2[1]}"
    return f"{b1}-{b2}"


def write_junction_reads_fasta(junction_reads, output_dir, bam_name):
    pair_to_fh = {}
    for r in junction_reads:
        key = pair_label(*r["breakpoints"])
        if key not in pair_to_fh:
            pair_to_fh[key] = open(os.path.join(output_dir, f"{bam_name}_{key}_junction_reads.fasta"), "w")
        pair_to_fh[key].write(f">{bam_name}:{r['read_name']}:{key}\n{r['capitalised_seq']}\n")
    for fh in pair_to_fh.values():
        fh.close()


def write_junction_reads_fastq(bam_file, junction_reads, output_dir, bam_name):
    """Write junction reads as FASTQ from the sequences kept during the search (*bam_file* is not reread)."""
    pair_to_fh = {}
    for r in junction_reads:
        key = pair_label(*r["breakpoints"])
        if key not in pair_to_fh:
            pair_to_fh[key] = open(os.path.join(output_dir, f"{bam_name}_{key}_junction_reads.fastq"), "w")
        seq = r["sequence"]
        quals = r["qualities"] if r["qualities"] is not None else [30] * len(seq)
        qual_str = "".join(chr(q + 33) for q in quals)
        pair_to_fh[key].write(f"@{bam_name}:{r['read_name']}:{key}\n{seq}\n+\n{qual_str}\n")
    for fh in pair_to_fh.values():
        fh.close()
