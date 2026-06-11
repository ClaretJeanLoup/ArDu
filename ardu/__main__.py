"""
Entry point for ``python -m ardu`` and the ``ardu`` console script.

Usage
-----
    ardu coverage  -b bams.txt -r regions.bed -n REF -o run1 [...]
    ardu softclips -b bams.txt -i run1_breakpoints.tsv -o run1 [...]
    ardu junctions -b bams.txt -i run1_breakpoints.tsv -o run1 [...]
    ardu parse     -i run1_coverage.tsv -o run1_rcn.tsv --RCN
"""

import argparse

from ardu.cli import coverage, softclips, junctions, parser


def main() -> None:
    ap = argparse.ArgumentParser(
        prog="ardu",
        description="ArDu — Architecture of Duplications toolkit",
    )
    sub = ap.add_subparsers(dest="command", required=True)

    coverage.register(sub)
    softclips.register(sub)
    junctions.register(sub)
    parser.register(sub)

    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
