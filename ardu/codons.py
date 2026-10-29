"""
Codon genotyping.

Covers:
    - the mutation file (one codon per line) and its mutation names (G119S, L1014S/F)
    - the codon of each read (or read pair) over a codon, translated
    - the number of copies of each amino acid, from read shares and the RCN
"""

from __future__ import annotations

import re
from typing import NamedTuple

import pysam
from Bio.Data.CodonTable import TranslationError
from Bio.Seq import Seq

from ardu import bam

_NAME = re.compile(r"^([A-Z*])(\d+)([A-Z*](?:/[A-Z*])*)$")
_STRANDS = {"forward": "forward", "+": "forward", "reverse": "reverse", "-": "reverse"}
CLASSES = ("S", "R", "O")


class Mutation(NamedTuple):
    chrom: str
    start: int              # 1-based first base of the codon on the reference
    stop: int               # start + 2
    gene: str               # locus name, as in the region file of ardu coverage
    name: str               # e.g. G119S or L1014S/F
    strand: str             # forward or reverse
    ref_aa: str
    alt_aas: frozenset[str]

    @property
    def id(self) -> str:
        return f"{self.name}:{self.gene}:{self.chrom}:{self.start}-{self.stop}:{self.strand}"


def parse_mutation_name(name: str) -> tuple[str, frozenset[str]]:
    """``G119S`` gives ``("G", {"S"})``; ``L1014S/F`` gives ``("L", {"S", "F"})``."""
    match = _NAME.match(name)
    if not match:
        raise ValueError(f"Mutation name '{name}' is not of the form G119S or L1014S/F")
    return match.group(1), frozenset(match.group(3).split("/"))


def read_mutations(path: str) -> list[Mutation]:
    """Read the tab-separated mutation file.

    Columns: chromosome, start, stop (1-based, inclusive, ``stop = start + 2``),
    gene, mutation name, strand (forward/reverse or +/-). Empty lines and lines
    starting with ``#`` are skipped.
    """
    mutations = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            if len(parts) < 6:
                raise ValueError(f"{path} line {n}: 6 tab-separated columns expected, got {len(parts)}")
            chrom, start, stop, gene, name, strand = parts[:6]
            if strand.lower() not in _STRANDS:
                raise ValueError(f"{path} line {n}: strand '{strand}' is not forward, reverse, + or -")
            start, stop = int(start), int(stop)
            if stop != start + 2:
                raise ValueError(f"{path} line {n}: a codon spans 3 bases (stop = start + 2), "
                                 f"got {start}-{stop}")
            ref_aa, alt_aas = parse_mutation_name(name)
            mutations.append(Mutation(chrom, start, stop, gene, name, _STRANDS[strand.lower()],
                                      ref_aa, alt_aas))
    return mutations


def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def translate(codon: str, table: int = 1) -> str:
    """Amino acid of *codon* in NCBI translation *table*; ``X`` when ambiguous."""
    try:
        return str(Seq(codon.upper()).translate(table=table))
    except TranslationError:
        return "X"


def reference_codon(fasta: pysam.FastaFile, mutation: Mutation) -> str:
    """Codon of the reference, read on the coding strand."""
    seq = fasta.fetch(mutation.chrom, mutation.start - 1, mutation.stop).upper()
    return reverse_complement(seq) if mutation.strand == "reverse" else seq


def read_codons(samfile: pysam.AlignmentFile, mutation: Mutation) -> tuple[dict[str, int], int]:
    """Count the codons read over *mutation*, on the coding strand.

    A codon is kept when its three bases pass the base quality filter, from one
    read or from the two mates of a pair (counted once). Where mates overlap,
    htslib gives the base to one mate only (quality summed when they agree, the
    better base down-weighted to 80% when they differ), so the base quality
    filter keeps at most one of them. Reads with a deletion or an insertion
    inside the codon are left out and returned as the second value.
    """
    positions = range(mutation.start - 1, mutation.start + 2)
    bases: dict[str, dict[int, tuple[str, int]]] = {}
    with_indel: set[str] = set()
    for column in samfile.pileup(
        mutation.chrom, positions[0], positions[-1] + 1, truncate=True, stepper="nofilter",
        min_base_quality=0, max_depth=10_000_000, ignore_overlaps=True,
    ):
        pos = column.reference_pos
        for pr in column.pileups:
            read = pr.alignment
            if read.is_supplementary or not bam.passes_read_filters(read):
                continue
            if pr.is_del or pr.is_refskip or (pr.indel > 0 and pos < positions[-1]):
                with_indel.add(read.query_name)
                continue
            quals = read.query_qualities
            qual = quals[pr.query_position] if quals is not None else bam.MIN_BASE_QUALITY
            if qual < bam.MIN_BASE_QUALITY:
                continue
            seen = bases.setdefault(read.query_name, {})
            if pos not in seen or qual > seen[pos][1]:
                seen[pos] = (read.query_sequence[pr.query_position].upper(), qual)

    counts: dict[str, int] = {}
    for name, seen in bases.items():
        if name in with_indel or len(seen) < 3:
            continue
        codon = "".join(seen[p][0] for p in positions)
        if mutation.strand == "reverse":
            codon = reverse_complement(codon)
        counts[codon] = counts.get(codon, 0) + 1
    return counts, len(with_indel)


def aa_class(aa: str, mutation: Mutation) -> str:
    """``S`` for the reference amino acid, ``R`` for a named mutant one, ``O`` otherwise."""
    if aa == mutation.ref_aa:
        return "S"
    return "R" if aa in mutation.alt_aas else "O"


def total_copies(rcn: float | None, locus_ploidy: int, norm_ploidy: int, dup_rcn: float) -> int:
    """Number of copies of the locus in a sample.

    The locus has ``rcn * norm_ploidy`` copies; it counts as duplicated when
    that is above ``dup_rcn`` times its ploidy, and otherwise keeps its ploidy,
    so that an RCN of 1.3 in a diploid gives 2 copies, not 3. Without an RCN
    the locus is taken as not duplicated.
    """
    if rcn is None:
        return locus_ploidy
    copies = rcn * norm_ploidy
    if copies <= dup_rcn * locus_ploidy:
        return locus_ploidy
    return max(locus_ploidy, round(copies))


def allocate_copies(
    reads: dict[str, int],
    n_copies: int,
    min_reads: int = 3,
    min_frac: float = 0.05,
) -> dict[str, int]:
    """Share *n_copies* between the amino acids of *reads* (``{aa: reads}``).

    An amino acid is kept when it has at least *min_reads* reads and at least
    half a copy's share of them, ``0.5 / n_copies``, or *min_frac* if higher.
    The copies are shared in proportion to the reads of the kept amino acids,
    by largest remainder, so they always sum to *n_copies*. A kept amino acid
    can get 0 copies when there are more of them than copies.
    """
    depth = sum(reads.values())
    if depth == 0 or n_copies < 1:
        return {}
    floor = max(min_frac, 0.5 / n_copies)
    kept = {aa: r for aa, r in reads.items() if r >= min_reads and r / depth >= floor}
    if not kept:
        return {}
    total = sum(kept.values())
    exact = {aa: r / total * n_copies for aa, r in kept.items()}
    copies = {aa: int(x) for aa, x in exact.items()}
    left = n_copies - sum(copies.values())
    # ties go to the amino acid with more reads, then alphabetically, so the result is stable
    for aa in sorted(exact, key=lambda a: (copies[a] - exact[a], -kept[a], a))[:left]:
        copies[aa] += 1
    return copies


def genotype_string(copies_by_class: dict[str, int]) -> str:
    """``{"S": 2, "R": 1}`` gives ``2S1R``."""
    return "".join(f"{copies_by_class[c]}{c}" for c in CLASSES if copies_by_class.get(c))
