"""
``ardu parse`` — extract individual metric columns from ArDu coverage output.

Corresponds to ArDuParser.py, refactored as a subcommand.
"""

from __future__ import annotations

import sys
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
                        help="Output TSV file (or header list if --format=list).")
    parser.add_argument("--uDOC",     action="store_true", help="Extract mean depth of coverage.")
    parser.add_argument("--sdDOC",    action="store_true", help="Extract standard deviation.")
    parser.add_argument("--medDOC",   action="store_true", help="Extract median depth of coverage.")
    parser.add_argument("--CovBases", action="store_true", help="Extract covered bases.")
    parser.add_argument("--RCN",      action="store_true", help="Extract relative copy number.")

    # --- filtering -----------------------------------------------------------
    parser.add_argument(
        "--filter",
        metavar="LOCUS:OP:VALUE",
        help=(
            "Filter columns by the value in a specific locus row. "
            "Format: LOCUS:OPERATOR:VALUE, e.g. 'chr1:100-200:>:1.4' or 'gene1:>=:2.0'. "
            "Supported operators: >, <, >=, <=, ==, !=. "
            "The locus column is always kept regardless of the filter result."
        ),
    )
    parser.add_argument(
        "--filter-locus",
        metavar="LOCUS",
        help="Locus row to filter on (alternative to embedding it in --filter).",
    )
    parser.add_argument(
        "--filter-op",
        metavar="OP",
        choices=[">", "<", ">=", "<=", "==", "!="],
        help="Comparison operator for --filter-locus / --filter-value.",
    )
    parser.add_argument(
        "--filter-value",
        metavar="VALUE",
        type=float,
        help="Numeric threshold for --filter-locus / --filter-op.",
    )

    # --- output format -------------------------------------------------------
    parser.add_argument(
        "--format",
        choices=["tsv", "list"],
        default="tsv",
        help=(
            "Output format: 'tsv' (default) writes the full filtered matrix; "
            "'list' writes one sample column header per line (useful for bam lists)."
        ),
    )
    parser.add_argument(
        "--suffix",
        default="",
        metavar="SUFFIX",
        help="Append a suffix to every header in list mode, e.g. '.bam'.",
    )

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_value(cell: str, field_index: int) -> str:
    """Extract a field from a semicolon-delimited coverage cell."""
    try:
        return str(cell).split(";")[field_index]
    except Exception:
        return "NA"


def _resolve_filter(args):
    """
    Return (locus, operator, value) from either --filter or the three
    individual --filter-locus / --filter-op / --filter-value flags.
    Returns (None, None, None) when no filter is requested.
    """
    if args.filter:
        # Expect "LOCUS:OP:VALUE" — locus itself may contain colons, so split
        # from the right to grab op and value first.
        parts = args.filter.rsplit(":", 2)
        if len(parts) != 3:
            raise ValueError(
                "--filter must be in the form LOCUS:OPERATOR:VALUE, "
                f"e.g. 'gene1:>:1.4'. Got: {args.filter!r}"
            )
        locus, op, raw_value = parts
        return locus, op, float(raw_value)

    if args.filter_locus or args.filter_op or args.filter_value is not None:
        missing = [
            flag for flag, val in (
                ("--filter-locus", args.filter_locus),
                ("--filter-op",    args.filter_op),
                ("--filter-value", args.filter_value),
            ) if val is None
        ]
        if missing:
            raise ValueError(
                f"All three flags are required together: {', '.join(missing)}"
            )
        return args.filter_locus, args.filter_op, float(args.filter_value)

    return None, None, None


_OPS = {
    ">":  lambda a, b: a > b,
    "<":  lambda a, b: a < b,
    ">=": lambda a, b: a >= b,
    "<=": lambda a, b: a <= b,
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
}


def _filter_columns(df: pd.DataFrame, locus: str, op: str, value: float) -> pd.DataFrame:
    """
    Keep only sample columns where the numeric value in *locus* row passes
    the comparison. The 'locus' column is always retained.
    """
    row = df[df["locus"] == locus]
    if row.empty:
        raise ValueError(f"Locus {locus!r} not found in the locus column.")

    compare = _OPS[op]
    keep = ["locus"]
    for col in df.columns:
        if col == "locus":
            continue
        try:
            cell_val = float(row[col].values[0])
        except (ValueError, TypeError):
            continue  # skip non-numeric / NA cells
        if compare(cell_val, value):
            keep.append(col)

    return df[keep]


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

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

    # Build the expanded metric matrix
    output_cols: dict = {"locus": df["locus"]}
    for sample in sample_columns:
        for field in requested_fields:
            idx = FIELD_MAP[field]
            output_cols[f"{sample}_{field}"] = df[sample].apply(
                lambda x, i=idx: _parse_value(x, i)
            )

    result = pd.DataFrame(output_cols)

    # Apply column filter if requested
    locus, op, value = _resolve_filter(args)
    if locus is not None:
        result = _filter_columns(result, locus, op, value)

    # Write output
    if args.format == "list":
        sample_headers = [c for c in result.columns if c != "locus"]
        lines = [f"{h}{args.suffix}" for h in sample_headers]
        with open(args.output, "w") as fh:
            fh.write("\n".join(lines) + "\n")
        print(f"Header list ({len(lines)} entries) saved to {args.output}")
    else:
        result.to_csv(args.output, sep="\t", index=False)
        print(f"Parsed output saved to {args.output}")