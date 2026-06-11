"""
``ardu parse`` — extract individual metric columns from ArDu coverage output.

Corresponds to ArDuParser.py, refactored as a subcommand.
"""

from __future__ import annotations

import pandas as pd


# Mapping of metric name → position in the semicolon-separated cell value
FIELD_MAP = {
    "uDOC":     0,   # mean depth of coverage
    "sdDOC":    1,   # standard deviation
    "medDOC":   2,   # median depth of coverage
    "CovBases": 3,   # covered bases
    "RCN":      4,   # relative copy number
}


# ---------------------------------------------------------------------------
# Argparse registration
# ---------------------------------------------------------------------------

def register(subparsers) -> None:
    """Add the ``parse`` sub-command to *subparsers*."""
    parser = subparsers.add_parser(
        "parse",
        help="Extract individual metric columns from an ArDu coverage TSV.",
        description=(
            "ArDu parse — convert the semicolon-packed ArDu coverage TSV "
            "into a flat, human-readable matrix of one or more metrics."
        ),
    )

    parser.add_argument("-i", "--input", required=True,
                        help="Input ArDu coverage TSV file.")
    parser.add_argument("-o", "--output", required=True,
                        help="Output TSV file.")
    parser.add_argument("--uDOC",     action="store_true", help="Extract mean depth of coverage.")
    parser.add_argument("--sdDOC",    action="store_true", help="Extract standard deviation.")
    parser.add_argument("--medDOC",   action="store_true", help="Extract median depth of coverage.")
    parser.add_argument("--CovBases", action="store_true", help="Extract covered bases.")
    parser.add_argument("--RCN",      action="store_true", help="Extract relative copy number.")

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def _parse_value(cell: str, field_index: int) -> str:
    """Extract a field from a semicolon-delimited coverage cell."""
    try:
        return str(cell).split(";")[field_index]
    except Exception:
        return "NA"


def run(args) -> None:
    """Entry point for ``ardu parse``."""
    requested_fields = [f for f in FIELD_MAP if getattr(args, f, False)]

    if not requested_fields:
        raise ValueError(
            "No fields selected. Use at least one of: "
            "--uDOC --sdDOC --medDOC --CovBases --RCN"
        )

    df = pd.read_csv(args.input, sep="\t", comment="#")
    sample_columns = [c for c in df.columns if c != "locus"]

    output_cols: dict = {"locus": df["locus"]}
    for sample in sample_columns:
        for field in requested_fields:
            idx = FIELD_MAP[field]
            output_cols[f"{sample}_{field}"] = df[sample].apply(
                lambda x, i=idx: _parse_value(x, i)
            )

    pd.DataFrame(output_cols).to_csv(args.output, sep="\t", index=False)
    print(f"Parsed output saved to {args.output}")
