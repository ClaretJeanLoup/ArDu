"""
Locus coordinates for the region file, from an annotation or from blast hits.

Covers:
    - reading a gff3 or gtf annotation and finding genes or transcripts by name
    - collecting the gene span, exons or CDS of a match, merged across isoforms
    - chaining the blast hits of a query into its best collinear set

All coordinates are 1-based and inclusive, as in gff, blast and the region file.
"""

from __future__ import annotations

import gzip
import math
import re
from dataclasses import dataclass, field
from urllib.parse import unquote

# Feature types that are never searched by name: parts of a transcript, and
# whole-sequence records that RefSeq gffs name after the chromosome.
LEAF_TYPES = {
    "exon", "cds", "five_prime_utr", "three_prime_utr", "utr", "5utr", "3utr",
    "start_codon", "stop_codon", "intron", "selenocysteine", "noncoding_exon",
}
SEQUENCE_TYPES = {"region", "chromosome", "contig", "supercontig", "scaffold", "sequence_feature"}

# Attributes searched for a name. Dbxref values are also searched without their
# database prefix (VectorBase:AGAP001356 matches AGAP001356).
NAME_KEYS = (
    "ID", "Name", "Alias", "gene", "gene_name", "gene_id", "gene_symbol", "locus_tag",
    "locus", "sequence_name", "transcript_id", "transcript_name", "Dbxref", "db_xref",
)

REGION_RE = re.compile(r"^(.+):([\d,]+)-([\d,]+)$")
VERSION_RE = re.compile(r"\.\d+$")


@dataclass
class Feature:
    chrom: str
    ftype: str
    start: int
    end: int
    strand: str
    fid: str
    parents: list[str] = field(default_factory=list)
    attrs: dict[str, list[str]] = field(default_factory=dict)

    def label(self) -> str:
        return f"{self.fid} ({self.ftype} {self.chrom}:{self.start}-{self.end})"


# ---------------------------------------------------------------------------
# Interval helpers
# ---------------------------------------------------------------------------

def merge_intervals(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Merge overlapping or touching 1-based intervals, sorted by chromosome and start."""
    merged: list[list] = []
    for chrom, start, end in sorted(intervals):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2] + 1:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([chrom, start, end])
    return [tuple(m) for m in merged]


def parse_region_line(text: str) -> tuple[str, int, int] | None:
    """``chrom:start-end`` (1-based) as ``(chrom, start, end)``, or None if *text* is not a region."""
    m = REGION_RE.match(text)
    if not m:
        return None
    start, end = int(m.group(2).replace(",", "")), int(m.group(3).replace(",", ""))
    if start < 1 or end < start:
        raise ValueError(f"Invalid region {text}: start must be >= 1 and <= end.")
    return m.group(1), start, end


# ---------------------------------------------------------------------------
# Annotation parsing
# ---------------------------------------------------------------------------

def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def _gff3_attrs(text: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = {}
    for part in text.strip().strip(";").split(";"):
        if "=" not in part:
            continue
        key, _, value = part.partition("=")
        attrs[key.strip()] = [unquote(v) for v in value.split(",")]
    return attrs


def _gtf_attrs(text: str) -> dict[str, list[str]]:
    attrs: dict[str, list[str]] = {}
    for part in text.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, value = part.partition(" ")
        attrs.setdefault(key, []).append(value.strip().strip('"'))
    return attrs


def _is_gtf(attr_text: str) -> bool:
    return re.match(r'^\s*\S+\s+"', attr_text) is not None


class Annotation:
    """A gff3 or gtf file held as features linked by their parents.

    In a gtf, genes and transcripts are built from ``gene_id`` and
    ``transcript_id`` when the file has no gene or transcript lines.
    """

    def __init__(self, path: str):
        self.features: dict[str, Feature] = {}
        self.children: dict[str, list[str]] = {}
        self.fmt = None
        self._read(path)
        for fid, feat in self.features.items():
            for p in feat.parents:
                if p in self.features:
                    self.children.setdefault(p, []).append(fid)
        self._index_names()

    def _read(self, path: str) -> None:
        with _open(path) as fh:
            for n, line in enumerate(fh):
                if line.startswith("##FASTA"):
                    break
                if not line.strip() or line.startswith("#"):
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) != 9:
                    raise ValueError(f"{path} line {n + 1}: expected 9 tab-separated columns.")
                if self.fmt is None:
                    self.fmt = "gtf" if _is_gtf(cols[8]) else "gff3"
                chrom, _, ftype, start, end, _, strand, _, attr_text = cols
                start, end = int(start), int(end)
                if self.fmt == "gtf":
                    self._add_gtf(chrom, ftype, start, end, strand, _gtf_attrs(attr_text), n)
                else:
                    attrs = _gff3_attrs(attr_text)
                    fid = attrs.get("ID", [f"_line{n + 1}"])[0]
                    if fid in self.features:        # one CDS on several lines shares its ID
                        fid = f"{fid}_line{n + 1}"
                    self.features[fid] = Feature(chrom, ftype, start, end, strand, fid,
                                                 attrs.get("Parent", []), attrs)

    def _add_gtf(self, chrom, ftype, start, end, strand, attrs, n) -> None:
        gene_id = attrs.get("gene_id", [None])[0]
        tx_id = attrs.get("transcript_id", [None])[0]
        gene_key = f"gene:{gene_id}" if gene_id else None
        tx_key = f"transcript:{tx_id}" if tx_id else None

        def node(key, kind, parents, keep):
            f = self.features.get(key)
            if f is None:
                kept = {k: v for k, v in attrs.items() if k in keep}
                f = self.features[key] = Feature(chrom, kind, start, end, strand, key, parents, kept)
            f.start, f.end = min(f.start, start), max(f.end, end)
            return f

        lower = ftype.lower()
        if lower == "gene" and gene_key:
            node(gene_key, "gene", [], attrs).attrs.update(attrs)
        elif lower in ("transcript", "mrna") and tx_key:
            node(tx_key, "transcript", [gene_key] if gene_key else [], attrs).attrs.update(attrs)
        else:
            if gene_key:
                node(gene_key, "gene", [], ("gene_id", "gene_name", "gene_symbol"))
            if tx_key:
                node(tx_key, "transcript", [gene_key] if gene_key else [],
                     ("gene_id", "gene_name", "transcript_id", "transcript_name"))
            parent = tx_key or gene_key
            fid = f"_line{n + 1}"
            self.features[fid] = Feature(chrom, ftype, start, end, strand, fid,
                                         [parent] if parent else [], attrs)

    def _index_names(self) -> None:
        self.names: dict[str, set[str]] = {}
        self.unversioned: dict[str, set[str]] = {}
        for fid, f in self.features.items():
            if f.ftype.lower() in LEAF_TYPES or f.ftype.lower() in SEQUENCE_TYPES:
                continue
            values = set()
            for key in NAME_KEYS:
                for v in f.attrs.get(key, []):
                    values.add(v)
                    if key in ("Dbxref", "db_xref") and ":" in v:
                        values.add(v.split(":", 1)[1])
            if f.fid.startswith(("gene:", "transcript:")):
                values.add(f.fid.split(":", 1)[1])
            for v in values:
                self.names.setdefault(v.lower(), set()).add(fid)
                unversioned = VERSION_RE.sub("", v)
                if unversioned != v:
                    self.unversioned.setdefault(unversioned.lower(), set()).add(fid)

    def _ancestors(self, fid: str) -> set[str]:
        seen: set[str] = set()
        todo = list(self.features[fid].parents)
        while todo:
            p = todo.pop()
            if p in seen or p not in self.features:
                continue
            seen.add(p)
            todo.extend(self.features[p].parents)
        return seen

    def find(self, name: str) -> list[Feature]:
        """Genes or transcripts matching *name* (any case), keeping only the outermost
        of nested matches (a gene, not also its mRNA carrying the same gene name).

        Without an exact match, IDs are also matched without their version or isoform
        number (F36A4.7 matches the WormBase transcripts F36A4.7.1 and F36A4.7.2).
        Isoforms of one gene that all match give that gene.
        """
        hits = self.names.get(name.lower()) or self.unversioned.get(name.lower(), set())
        top = {h for h in hits if not (self._ancestors(h) & hits)}
        parents = {p for h in top for p in self.features[h].parents}
        if len(top) > 1 and len(parents) == 1 and all(self.features[h].parents for h in top):
            parent = parents.pop()
            if parent in self.features:
                top = {parent}
        return sorted((self.features[h] for h in top), key=lambda f: (f.chrom, f.start))

    def intervals(self, feat: Feature, feature: str) -> list[tuple[str, int, int]]:
        """Span of *feat* (``feature='gene'``) or its exons or CDS, merged across isoforms."""
        if feature == "gene":
            return [(feat.chrom, feat.start, feat.end)]
        out, todo, seen = [], [feat.fid], set()
        while todo:
            fid = todo.pop()
            if fid in seen:
                continue
            seen.add(fid)
            f = self.features[fid]
            if f.ftype.lower() == feature.lower():
                out.append((f.chrom, f.start, f.end))
            todo.extend(self.children.get(fid, []))
        return merge_intervals(out)


# ---------------------------------------------------------------------------
# Blast hit chaining
# ---------------------------------------------------------------------------

@dataclass
class Hsp:
    chrom: str
    strand: str
    qstart: int
    qend: int
    sstart: int          # sstart <= send whatever the strand
    send: int
    bits: float
    pident: float
    length: int


@dataclass
class Chain:
    hsps: list[Hsp]
    score: float

    @property
    def chrom(self) -> str:
        return self.hsps[0].chrom

    @property
    def strand(self) -> str:
        return self.hsps[0].strand

    @property
    def start(self) -> int:
        return min(h.sstart for h in self.hsps)

    @property
    def end(self) -> int:
        return max(h.send for h in self.hsps)

    def query_cover(self, qlen: int) -> float:
        covered = merge_intervals([("q", h.qstart, h.qend) for h in self.hsps])
        return sum(e - s + 1 for _, s, e in covered) / qlen

    def pident(self) -> float:
        total = sum(h.length for h in self.hsps)
        return sum(h.pident * h.length for h in self.hsps) / total


def hsps_from_blast(fields: list[str], min_intron: int) -> list[Hsp]:
    """Hits from a blast tabular line with the columns ``qseqid sseqid pident length
    qstart qend sstart send evalue bitscore qlen btop``.

    A gapped hit can bridge a short intron (C. elegans has many of 45 to 50 bp), so the
    hit is cut wherever the reference has a run of at least *min_intron* bases missing
    from the query; the bitscore is shared by query length.
    """
    qstart, qend, sstart, send = (int(x) for x in fields[4:8])
    bits, pident = float(fields[9]), float(fields[2])
    strand, step = ("+", 1) if sstart <= send else ("-", -1)

    pieces = []                   # (qstart, qend, sstart, send) in alignment order
    q, s = qstart, sstart
    cq, cs, last_q, last_s, run = qstart, sstart, None, None, 0
    for token in re.findall(r"\d+|\D\D", fields[11]):
        columns = [("M", "M")] * int(token) if token.isdigit() else [(token[0], token[1])]
        for qc, sc in columns:
            if qc == "-":                     # base of the reference only
                run += 1
                s += step
                continue
            if run >= min_intron and last_q is not None:
                pieces.append((cq, last_q, cs, last_s))
                cq, cs = q, s
            run = 0
            last_q = q
            q += 1
            if sc != "-":
                last_s = s
                s += step
    pieces.append((cq, last_q, cs, last_s))

    total = qend - qstart + 1
    return [Hsp(fields[1], strand, a, b, min(c, d), max(c, d),
                bits * (b - a + 1) / total, pident, b - a + 1)
            for a, b, c, d in pieces]


# Chain score costs. Blast extends hits a few bases past exon edges by chance, so the
# query bases two hits share are charged about what they earned (megablast gives about
# 1.85 bits per matching base); otherwise these chance extensions would decide between
# two copies of a gene. The gap cost only breaks ties in favour of the compact chain.
OVERLAP_BITS = 2.0
GAP_BITS = 0.1          # per doubling of the gap


def _link_cost(prev: Hsp, nxt: Hsp, max_gap: int, overlap: int) -> float | None:
    """Cost of putting *nxt* after *prev* in a chain, or None when it cannot follow:
    same chromosome and strand, further along the query, and further along the genome
    in the query's orientation, at most *max_gap* bp away. Overlaps of up to *overlap*
    bp are allowed, on the query and on the genome."""
    if nxt.chrom != prev.chrom or nxt.strand != prev.strand:
        return None
    q_overlap = prev.qend - nxt.qstart + 1
    if nxt.qstart <= prev.qstart or nxt.qend <= prev.qend or q_overlap > overlap:
        return None
    if prev.strand == "+":
        if nxt.sstart <= prev.sstart or nxt.send <= prev.send:
            return None
        gap = nxt.sstart - prev.send - 1
    else:
        if nxt.send >= prev.send or nxt.sstart >= prev.sstart:
            return None
        gap = prev.sstart - nxt.send - 1
    if not -overlap <= gap <= max_gap:
        return None
    return OVERLAP_BITS * max(0, q_overlap) + GAP_BITS * math.log2(max(gap, 0) + 2)


def best_chain(hsps: list[Hsp], max_gap: int, overlap: int) -> Chain | None:
    """The collinear set of hits with the highest summed bitscore, less the link costs."""
    if not hsps:
        return None
    hsps = sorted(hsps, key=lambda h: (h.qstart, h.qend))
    score = [h.bits for h in hsps]
    back = [-1] * len(hsps)
    for i, h in enumerate(hsps):
        for j in range(i):
            cost = _link_cost(hsps[j], h, max_gap, overlap)
            if cost is not None and score[j] + h.bits - cost > score[i]:
                score[i], back[i] = score[j] + h.bits - cost, j
    i = max(range(len(hsps)), key=score.__getitem__)
    best = score[i]
    chain = []
    while i != -1:
        chain.append(hsps[i])
        i = back[i]
    return Chain(chain[::-1], best)


def best_two_chains(hsps: list[Hsp], max_gap: int, overlap: int,
                    max_hsps: int = 500) -> tuple[Chain | None, Chain | None]:
    """The best chain, and the best chain among the hits outside its span (a second
    copy of the locus in the reference, or a paralogue)."""
    hsps = sorted(hsps, key=lambda h: -h.bits)[:max_hsps]
    first = best_chain(hsps, max_gap, overlap)
    if first is None:
        return None, None
    rest = [h for h in hsps if not (h.chrom == first.chrom
                                    and h.sstart <= first.end and h.send >= first.start)]
    return first, best_chain(rest, max_gap, overlap)


def alternative_hits(chain: Chain, hsps: list[Hsp], ratio: float) -> list[tuple[Hsp, Hsp]]:
    """``(chain hit, other hit)`` pairs where a hit outside *chain* covers at least half of
    the same query bases, elsewhere on the genome, with at least *ratio* of its bitscore.
    This catches a second copy close enough to the first that one chain can mix them."""
    pairs = []
    for h in chain.hsps:
        for g in hsps:
            if g in chain.hsps or (g.chrom == h.chrom and g.sstart <= h.send and g.send >= h.sstart):
                continue
            shared = min(h.qend, g.qend) - max(h.qstart, g.qstart) + 1
            if shared >= 0.5 * (h.qend - h.qstart + 1) and g.bits >= ratio * h.bits:
                pairs.append((h, g))
                break
    return pairs
