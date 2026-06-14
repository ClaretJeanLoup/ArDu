"""
``ardu parse`` — extract individual metric columns from ArDu coverage output.
"""

from __future__ import annotations

import os
import pandas as pd

FIELD_MAP = {
    "uDOC":     0,
    "sdDOC":    1,
    "medDOC":   2,
    "CovBases": 3,
    "RCN":      4,
}

_OPS = {
    ">":  lambda a, b: a > b,
    "<":  lambda a, b: a < b,
    ">=": lambda a, b: a >= b,
    "<=": lambda a, b: a <= b,
    "==": lambda a, b: a == b,
    "!=": lambda a, b: a != b,
}

_OP_LABELS = {">": "gt", "<": "lt", ">=": "gte", "<=": "lte", "==": "eq", "!=": "neq"}


# ---------------------------------------------------------------------------
# Argparse registration
# ---------------------------------------------------------------------------

def register(subparsers) -> None:
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
    parser.add_argument("-o", "--output", default=None, metavar="PREFIX",
                        help=(
                            "Optional output prefix. The file is always written next to "
                            "the input file. If omitted, the filename is built from the "
                            "selected fields and filter parameters. "
                            "e.g. -o mystudy → mystudy_RCN_ace1_gt_1.4.txt"
                        ))
    parser.add_argument("--uDOC",     action="store_true", help="Extract mean depth of coverage.")
    parser.add_argument("--sdDOC",    action="store_true", help="Extract standard deviation.")
    parser.add_argument("--medDOC",   action="store_true", help="Extract median depth of coverage.")
    parser.add_argument("--CovBases", action="store_true", help="Extract covered bases.")
    parser.add_argument("--RCN",      action="store_true", help="Extract relative copy number.")

    # filtering
    parser.add_argument(
        "--filter", metavar="LOCUS:OP:VALUE",
        help="Filter columns by value at a locus. Format: LOCUS:OP:VALUE e.g. 'ace1:>:1.4'.",
    )
    parser.add_argument("--filter-locus", metavar="LOCUS")
    parser.add_argument("--filter-op",    metavar="OP", choices=list(_OPS))
    parser.add_argument("--filter-value", metavar="VALUE", type=float)

    # output format
    parser.add_argument(
        "--format", choices=["tsv", "list"], default="tsv",
        help="'tsv' — full matrix (default); 'list' — one header per line.",
    )
    parser.add_argument(
        "--suffix", default=None, metavar="SUFFIX",
        help=(
            "Suffix appended to each entry in list mode. "
            "Defaults to '.bam' when --format=list."
        ),
    )

    parser.set_defaults(func=run)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_value(cell: str, field_index: int) -> str:
    try:
        return str(cell).split(";")[field_index]
    except Exception:
        return "NA"


def _resolve_filter(args):
    """Return (locus, operator, value) or (None, None, None)."""
    if args.filter:
        parts = args.filter.rsplit(":", 2)
        if len(parts) != 3:
            raise ValueError(
                f"--filter must be LOCUS:OP:VALUE. Got: {args.filter!r}"
            )
        locus, op, raw_value = parts
        return locus, op, float(raw_value)

    if any(v is not None for v in (args.filter_locus, args.filter_op, args.filter_value)):
        missing = [
            flag for flag, val in (
                ("--filter-locus", args.filter_locus),
                ("--filter-op",    args.filter_op),
                ("--filter-value", args.filter_value),
            ) if val is None
        ]
        if missing:
            raise ValueError(f"All three flags required together: {', '.join(missing)}")
        return args.filter_locus, args.filter_op, float(args.filter_value)

    return None, None, None


def _filter_columns(df: pd.DataFrame, locus: str, op: str, value: float) -> pd.DataFrame:
    row = df[df["locus"] == locus]
    if row.empty:
        raise ValueError(f"Locus {locus!r} not found.")
    compare = _OPS[op]
    keep = ["locus"]
    for col in df.columns:
        if col == "locus":
            continue
        try:
            if compare(float(row[col].values[0]), value):
                keep.append(col)
        except (ValueError, TypeError):
            pass
    return df[keep]


def _auto_output(input_path: str, prefix: str | None, fields: list[str],
                 locus, op, value, fmt: str) -> str:
    """
    Build output path:
    - Directory is always taken from input file location
    - Filename: [prefix_]fields[_locus_op_value].ext
    - ext is .txt for list mode, .tsv for tsv mode
    """
    input_dir = os.path.dirname(os.path.abspath(input_path))

    parts = []
    if prefix:
        parts.append(prefix)
    parts.extend(fields)
    if locus is not None:
        parts.append(f"{locus}_{_OP_LABELS[op]}_{value:g}")

    ext = ".txt" if fmt == "list" else ".tsv"
    return os.path.join(input_dir, "_".join(parts) + ext)


# ---------------------------------------------------------------------------
# Main orchestration
# ---------------------------------------------------------------------------

def run(args) -> None:
    requested_fields = [f for f in FIELD_MAP if getattr(args, f, False)]
    if not requested_fields:
        raise ValueError(
            "No fields selected. Use at least one of: "
            "--uDOC --sdDOC --medDOC --CovBases --RCN"
        )

    locus, op, value = _resolve_filter(args)

    output_path = _auto_output(
        args.input, args.output, requested_fields, locus, op, value, args.format
    )

    suffix = args.suffix if args.suffix is not None else (".bam" if args.format == "list" else "")

    df = pd.read_csv(args.input, sep="\t", comment="#")
    sample_columns = [c for c in df.columns if c != "locus"]

    output_cols: dict = {"locus": df["locus"]}
    for sample in sample_columns:
        for field in requested_fields:
            idx = FIELD_MAP[field]
            col_name = f"{sample}_{field}" if args.format == "tsv" else sample
            output_cols[col_name] = df[sample].apply(
                lambda x, i=idx: _parse_value(x, i)
            )

    result = pd.DataFrame(output_cols)

    if locus is not None:
        result = _filter_columns(result, locus, op, value)

    if args.format == "list":
        entries = [f"{c}{suffix}" for c in result.columns if c != "locus"]
        with open(output_path, "w") as fh:
            fh.write("\n".join(entries) + "\n")
        print(f"Header list ({len(entries)} entries) saved to {output_path}")
    else:
        result.to_csv(output_path, sep="\t", index=False)
        print(f"Parsed output saved to {output_path}")