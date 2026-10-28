"""Tests of ``ardu regions``: gff3 and gtf lookup, and the blast hit chains."""
from __future__ import annotations

import gzip
import random
import shutil
import subprocess
import sys

import pysam
import pytest

from ardu.intervals import get_regions
from ardu.locate import Annotation, Hsp, best_chain, best_two_chains, hsps_from_blast

# Gene g1 on chr1: two isoforms whose exons overlap, a CDS written over two lines with
# one ID (as in RefSeq), and a Dbxref. Gene g2 on chr2 shares the Name "dup" with g3.
GFF3 = """##gff-version 3
chr1\tsrc\tregion\t1\t5000\t.\t+\t.\tID=chr1;Name=g1
chr1\tsrc\tgene\t1000\t2000\t.\t+\t.\tID=gene-g1;Name=ABC1;Dbxref=VectorBase:AGAP1
chr1\tsrc\tmRNA\t1000\t2000\t.\t+\t.\tID=rna-g1a;Parent=gene-g1;gene=ABC1
chr1\tsrc\texon\t1000\t1100\t.\t+\t.\tID=e1;Parent=rna-g1a;gene=ABC1
chr1\tsrc\texon\t1500\t1600\t.\t+\t.\tID=e2;Parent=rna-g1a;gene=ABC1
chr1\tsrc\tCDS\t1050\t1100\t.\t+\t0\tID=cds-g1a;Parent=rna-g1a
chr1\tsrc\tCDS\t1500\t1550\t.\t+\t0\tID=cds-g1a;Parent=rna-g1a
chr1\tsrc\tmRNA\t1000\t2000\t.\t+\t.\tID=rna-g1b;Parent=gene-g1;gene=ABC1
chr1\tsrc\texon\t1000\t1120\t.\t+\t.\tID=e3;Parent=rna-g1b
chr1\tsrc\texon\t1900\t2000\t.\t+\t.\tID=e4;Parent=rna-g1b
chr2\tsrc\tgene\t100\t300\t.\t-\t.\tID=gene-g2;Name=dup
chr2\tsrc\tgene\t800\t900\t.\t-\t.\tID=gene-g3;Name=dup
"""

GTF = """chr1\tsrc\texon\t1000\t1100\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "abc1";
chr1\tsrc\texon\t1500\t1600\t.\t+\t.\tgene_id "G1"; transcript_id "T1"; gene_name "abc1";
chr1\tsrc\texon\t1000\t1120\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "abc1";
chr1\tsrc\texon\t1900\t2000\t.\t+\t.\tgene_id "G1"; transcript_id "T2"; gene_name "abc1";
"""


@pytest.fixture
def gff(tmp_path):
    path = tmp_path / "genes.gff3.gz"
    with gzip.open(path, "wt") as fh:
        fh.write(GFF3)
    return str(path)


def run_regions(cwd, *args):
    cmd = [sys.executable, "-m", "ardu", "regions", *args, "-o", "out"]
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def write_list(tmp_path, text):
    path = tmp_path / "loci.txt"
    path.write_text(text)
    return str(path)


def test_exons_are_merged_across_isoforms(gff):
    ann = Annotation(gff)
    (g1,) = ann.find("abc1")
    assert g1.fid == "gene-g1"
    assert ann.intervals(g1, "exon") == [("chr1", 1000, 1120), ("chr1", 1500, 1600),
                                         ("chr1", 1900, 2000)]
    assert ann.intervals(g1, "CDS") == [("chr1", 1050, 1100), ("chr1", 1500, 1550)]
    assert ann.intervals(g1, "gene") == [("chr1", 1000, 2000)]


def test_names_ids_transcripts_and_dbxref(gff):
    ann = Annotation(gff)
    assert [f.fid for f in ann.find("AGAP1")] == ["gene-g1"]
    assert [f.fid for f in ann.find("rna-g1b")] == ["rna-g1b"]
    assert ann.intervals(ann.find("rna-g1b")[0], "exon") == [("chr1", 1000, 1120),
                                                             ("chr1", 1900, 2000)]
    assert ann.find("g1") == []                              # the chr1 region record's Name
    assert len(ann.find("dup")) == 2


def test_gtf_without_gene_lines(tmp_path):
    path = tmp_path / "genes.gtf"
    path.write_text(GTF)
    ann = Annotation(str(path))
    assert ann.fmt == "gtf"
    (g,) = ann.find("ABC1")
    assert (g.start, g.end) == (1000, 2000)
    assert ann.intervals(g, "exon") == [("chr1", 1000, 1120), ("chr1", 1500, 1600),
                                        ("chr1", 1900, 2000)]
    assert ann.intervals(ann.find("T1")[0], "exon") == [("chr1", 1000, 1100),
                                                        ("chr1", 1500, 1600)]


def test_ids_without_version_or_isoform_number(tmp_path):
    path = tmp_path / "wb.gtf"
    path.write_text(
        'IV\tWormBase\texon\t100\t200\t.\t+\t.\tgene_id "WBGene1"; transcript_id "F36A4.7.1";\n'
        'IV\tWormBase\texon\t150\t300\t.\t+\t.\tgene_id "WBGene1"; transcript_id "F36A4.7.2";\n'
        'IV\tWormBase\texon\t900\t950\t.\t+\t.\tgene_id "WBGene2"; transcript_id "F36A4.71.1";\n')
    ann = Annotation(str(path))
    assert [f.fid for f in ann.find("F36A4.7")] == ["gene:WBGene1"]      # both isoforms
    assert [f.fid for f in ann.find("F36A4.7.2")] == ["transcript:F36A4.7.2"]
    assert [f.fid for f in ann.find("F36A4.71")] == ["transcript:F36A4.71.1"]   # one isoform
    assert ann.find("F36A4") == []


def test_region_file_pools_pads_and_maps(gff, tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">CHR1\n" + "A" * 2010 + "\n")
    pysam.faidx(str(fa))
    (tmp_path / "map.txt").write_text("chr1 CHR1\n")
    loci = write_list(tmp_path, "ABC1 locusA\nCHR1:1,150-1,300 locusA\n# comment\n"
                                "CHR1:10-20\n")
    out = run_regions(tmp_path, "-g", gff, "-l", loci, "--pad", "20", "--chrom-map",
                      "map.txt", "--fasta", str(fa), "--plot-margin", "500")
    assert out.returncode == 0, out.stderr
    regions = get_regions(str(tmp_path / "out_regions.bed"))
    # exons 1000-1120, 1500-1600 and 1900-2000 padded by 20, the region joins the first
    # exon after padding, the last one is clipped at the end of CHR1 (2010 bp).
    assert regions == {"locusA": ["CHR1:980-1320", "CHR1:1480-1620", "CHR1:1880-2010"],
                       "CHR1:10-20": ["CHR1:1-40"]}
    plot = (tmp_path / "out_plot.bed").read_text().splitlines()
    assert plot[0] == "locusA\tCHR1:480-2510"


def test_missing_and_ambiguous_names_stop(gff, tmp_path):
    loci = write_list(tmp_path, "ABC1\nnope\ndup\n")
    out = run_regions(tmp_path, "-g", gff, "-l", loci)
    assert out.returncode != 0
    assert "nope: not found" in out.stderr
    assert "dup: matches 2 features" in out.stderr and "gene-g3" in out.stderr
    assert not (tmp_path / "out_regions.bed").exists()


def test_unknown_chromosome_stops(gff, tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">other\nACGT\n")
    pysam.faidx(str(fa))
    out = run_regions(tmp_path, "-g", gff, "-l", write_list(tmp_path, "ABC1\n"),
                      "--fasta", str(fa))
    assert out.returncode != 0 and "--chrom-map" in out.stderr


# ---------------------------------------------------------------------------
# Hit chains
# ---------------------------------------------------------------------------

def hsp(chrom, strand, qs, qe, ss, se, bits):
    return Hsp(chrom, strand, qs, qe, ss, se, bits, 100.0, qe - qs + 1)


def test_chain_keeps_collinear_exons_over_a_stronger_decoy():
    exons = [hsp("c", "+", 1, 100, 1000, 1099, 180), hsp("c", "+", 101, 200, 3000, 3099, 180),
             hsp("c", "+", 201, 300, 5000, 5099, 180)]
    decoy = hsp("c", "+", 101, 200, 900, 999, 190)          # before exon 1 on the genome
    chain = best_chain(exons + [decoy], max_gap=10_000, overlap=30)
    assert chain.hsps == exons and chain.score == pytest.approx(540, abs=3)


def test_chain_on_the_minus_strand_and_gap_limit():
    hits = [hsp("c", "-", 1, 100, 5000, 5099, 180), hsp("c", "-", 101, 200, 3000, 3099, 180)]
    assert len(best_chain(hits, max_gap=10_000, overlap=30).hsps) == 2
    assert len(best_chain(hits, max_gap=1_000, overlap=30).hsps) == 1


def test_gapped_hit_is_cut_at_an_intron():
    # query 1-60 on the reference 101-205: 30 bp, a 45 bp intron, 30 bp with a
    # 3 bp deletion in the reference (kept in one piece).
    btop = "30" + "-A" * 45 + "10" + "C-" * 3 + "17"
    fields = ["q", "c", "99", "105", "1", "60", "101", "205", "1e-20", "100", "60", btop]
    plus = hsps_from_blast(fields, min_intron=20)
    assert [(h.qstart, h.qend, h.sstart, h.send) for h in plus] == [(1, 30, 101, 130),
                                                                    (31, 60, 176, 202)]
    assert sum(h.bits for h in plus) == pytest.approx(100)
    fields[6:8] = ["205", "101"]
    minus = hsps_from_blast(fields, min_intron=20)
    assert [(h.sstart, h.send, h.strand) for h in minus] == [(176, 205, "-"), (104, 130, "-")]
    assert len(hsps_from_blast(fields, min_intron=50)) == 1


def test_second_chain_is_another_copy():
    copy1 = [hsp("c", "+", 1, 100, 1000, 1099, 180), hsp("c", "+", 101, 200, 2000, 2099, 180)]
    copy2 = [hsp("d", "-", 1, 100, 9000, 9099, 170), hsp("d", "-", 101, 200, 8000, 8099, 170)]
    first, second = best_two_chains(copy1 + copy2, max_gap=10_000, overlap=30)
    assert first.hsps == copy1
    assert (second.chrom, second.start, second.end, second.strand) == ("d", 8000, 9099, "-")


# ---------------------------------------------------------------------------
# Blast end to end
# ---------------------------------------------------------------------------

@pytest.mark.skipif(shutil.which("blastn") is None, reason="blastn not installed")
def test_blast_finds_the_exons_and_flags_a_second_copy(tmp_path):
    rng = random.Random(3)
    seq = [rng.choice("ACGT") for _ in range(30_000)]
    ex1, ex2 = "".join(seq[2000:2300]), "".join(seq[4000:4400])      # 2,001-2,300, 4,001-4,400
    single = "".join(seq[10000:10500])                              # 10,001-10,500
    seq[20000:20300], seq[22000:22400] = list(ex1), list(ex2)       # second copy at 20,001
    (tmp_path / "ref.fa").write_text(">c\n" + "".join(seq) + "\n")
    (tmp_path / "q.fa").write_text(f">geneA\n{ex1}{ex2}\n>geneB\n{single}\n")

    out = run_regions(tmp_path, "-q", "q.fa", "--reference", "ref.fa")
    assert out.returncode == 0, out.stderr
    regions = get_regions(str(tmp_path / "out_regions.bed"))
    assert regions["geneB"] == ["c:10001-10500"]
    # blast may extend a hit a few bases past an exon edge by chance
    found = [tuple(map(int, r.split(":")[1].split("-"))) for r in regions["geneA"]]
    assert any(len(found) == 2 and all(abs(a - b) <= 10 for f, c in zip(found, copy)
                                       for a, b in zip(f, c))
               for copy in ([(2001, 2300), (4001, 4400)], [(20001, 20300), (22001, 22400)]))
    hits = {l.split("\t")[0]: l.rstrip("\n").split("\t")
            for l in open(tmp_path / "out_hits.tsv")}
    assert hits["geneA"][12] == "second_copy" and hits["geneB"][12] == "."
    assert "second chain" in out.stdout
