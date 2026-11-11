"""Tests of the overview depth profile: depth of each segment relative to the
flanks, numbered breakpoints, target gene and marked positions.

Four breakpoints cut the locus into flank 0, segments 1 to 3 and flank 4.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ardu.cli.junctions import _breakpoint_order, _depth_levels, _draw_profile, _read_marks, _read_target

BKPS = [1000, 2000, 3000, 4000]


def _segments(depths):
    bounds = [500] + BKPS + [4500]
    labels = ["flank", "1-2", "2-3", "3-4", "flank"]
    return [{"bam": "s", "locus": "gene", "segment": lab, "chrom": "2R", "start": a + 1, "end": b,
             "median_depth": d * 20, "relative_depth": d}
            for lab, a, b, d in zip(labels, bounds, bounds[1:], depths)]


def test_breakpoint_order_from_segments():
    assert _breakpoint_order(_segments([1] * 5)) == [1, 2, 3, 4]


def test_profile_figure_with_target_and_mark():
    segments = _segments([1, 1.9, 1.4, 1.9, 1])
    target = [("2R", 1200, 1300, "gene"), ("2R", 1400, 1600, "gene"), ("2R", 3200, 3300, "other")]
    fig, ax = plt.subplots()
    _draw_profile(ax, "s", "gene", segments, target, [("2R", 1450, "G119S")])
    texts = [t.get_text() for t in ax.texts]
    assert {"1.9", "1.4", "1.0"} <= set(texts)          # depth of each segment
    assert {"1", "2", "3", "4"} <= set(texts)              # numbered breakpoints
    assert "gene" in texts and "G119S" in texts
    # only the rows named after the locus are drawn, and no junction rows
    assert len(ax.patches) == 2
    assert not any("deletion" in x or "duplication" in x for x in texts)
    plt.close(fig)


def test_read_target_and_marks(tmp_path):
    bed = tmp_path / "loci.bed"
    bed.write_text("2R\t10\t20\tgene\n\n2R\t30\t40\n")
    mut = tmp_path / "mut.tsv"
    mut.write_text("2R\t15\tG119S\n2R\t35\n")
    assert _read_target(str(bed)) == [("2R", 10, 20, "gene"), ("2R", 30, 40, "")]
    assert _read_marks(str(mut)) == [("2R", 15, "G119S"), ("2R", 35, "")]


GFF = """##gff-version 3
NC_1\tRefSeq\tgene\t1100\t1700\t.\t+\t.\tID=g1;Name=ace-1
NC_1\tRefSeq\tmRNA\t1100\t1700\t.\t+\t.\tID=t1;Parent=g1
NC_1\tRefSeq\tCDS\t1200\t1300\t.\t+\t0\tID=c1;Parent=t1;product=acetylcholinesterase
NC_1\tRefSeq\tCDS\t1400\t1600\t.\t+\t0\tID=c1;Parent=t1;product=acetylcholinesterase
NC_1\tRefSeq\tgene\t2500\t3500\t.\t-\t.\tID=g2;locus_tag=AGAP9
NC_1\tRefSeq\tmRNA\t2500\t3500\t.\t-\t.\tID=t2;Parent=g2
NC_1\tRefSeq\tCDS\t2600\t2700\t.\t-\t0\tID=c2;Parent=t2;product=uncharacterised protein
NC_1\tRefSeq\tgene\t9000\t9500\t.\t+\t.\tID=g3;Name=far
"""


def test_gene_track_names_cds_and_window(tmp_path):
    from ardu.locate import GeneTrack
    gff = tmp_path / "a.gff"
    gff.write_text(GFF)
    genes = GeneTrack(str(gff), "CDS", "name", {"NC_1": "2R"}).genes("2R", 501, 4500, "ace-1")
    assert [(g["name"], g["target"], g["blocks"]) for g in genes] == [
        ("ace-1", True, [(1200, 1300), (1400, 1600)]), ("AGAP9", False, [(2600, 2700)])]
    product = GeneTrack(str(gff), "gene", "product", {"NC_1": "2R"}).genes("2R", 501, 4500)
    assert [g["name"] for g in product] == ["acetylcholinesterase", "uncharacterised protein"]
    assert product[0]["blocks"] == [(1100, 1700)]


def test_profile_figure_with_gff_genes(tmp_path):
    from ardu.locate import GeneTrack
    gff = tmp_path / "a.gff"
    gff.write_text(GFF)
    fig, ax = plt.subplots()
    _draw_profile(ax, "s", "ace-1", _segments([1, 1.9, 1.4, 1.9, 1]), None, None,
                  GeneTrack(str(gff), "CDS", "name", {"NC_1": "2R"}))
    texts = {t.get_text(): t for t in ax.texts}
    assert "ace-1" in texts and "far" not in texts
    assert texts["ace-1"].get_color() == "black"           # the target is labelled in black, by name
    # the other genes are numbered on the track and named in the legend
    assert "AGAP9" not in texts
    assert [t.get_text() for t in ax.get_legend().get_texts()] == ["1  AGAP9"]
    assert len(ax.patches) == 3                             # 2 CDS of ace-1, 1 of AGAP9
    plt.close(fig)


def test_small_depth_differences_are_one_level():
    # 1.1, 1.0 and 1.2 (within 0.2) are one level; 2.6 is a step
    x = [0, 1, 2, 3, 5]
    levels = _depth_levels([1.1, 1.0, 1.2, 2.6], x)
    assert [(a, b) for a, b, _ in levels] == [(0, 2), (3, 3)]
    assert abs(levels[0][2] - 1.1) < 1e-9 and levels[1][2] == 2.6
    assert _depth_levels([1.0, None, 1.0], [0, 1, 2, 3]) == [(0, 0, 1.0), (1, 1, None), (2, 2, 1.0)]


def test_mostly_uncovered_segment_is_dotted_and_explained():
    segments = _segments([1, 2.0, 1.0, 2.0, 1])
    segments[2]["pct_zero"] = 55.0
    fig, ax = plt.subplots()
    _draw_profile(ax, "s", "gene", segments)
    from ardu.cli.junctions import _GENE_COLOUR
    dotted = [l for l in ax.lines if l.get_color() == _GENE_COLOUR and len(l.get_xdata()) == 2
              and l.get_ydata()[0] == l.get_ydata()[1] and l.get_ydata()[0] > 0]
    assert len(dotted) == 1                                   # segment 2-3 only
    legend_texts = [t.get_text() for leg in ax.artists + [ax.get_legend()] if leg is not None
                    for t in leg.get_texts()]
    assert any("not reliable" in t for t in legend_texts)
    assert "scale is not reliable" not in ax.get_title()     # the flanks are covered
    plt.close(fig)
