#!/usr/bin/env python3

import argparse
import pandas as pd


# Mapping of fields to positions in the ArDu coverage outgile
FIELD_MAP = {
    "uDOC": 0,       # mean
    "sdDOC": 1,      # sd
    "medDOC": 2,     # median
    "CovBases": 3,   # coveredbases
    "RCN": 4         # normalised
}


def parse_value(cell, field_index):
    """Extract a specific value from the semicolon-separated coverage field."""
    try:
        parts = str(cell).split(";")
        return parts[field_index]
    except Exception:
        return "NA"


def main():

    parser = argparse.ArgumentParser(description="Parse ArDu coverage.tsv output into a simplified tab-delimited matrix.")

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input ArDu coverage.tsv file")

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output TSV file")

    parser.add_argument("--uDOC", action="store_true", help="Extract mean depth of coverage")
    parser.add_argument("--sdDOC", action="store_true", help="Extract standard deviation")
    parser.add_argument("--medDOC", action="store_true", help="Extract median depth of coverage")
    parser.add_argument("--CovBases", action="store_true", help="Extract covered bases")
    parser.add_argument("--RCN", action="store_true", help="Extract relative copy number")

    args = parser.parse_args()

    requested_fields = []

    for field in FIELD_MAP:
        if getattr(args, field):
            requested_fields.append(field)

    if not requested_fields:
        raise ValueError("No fields selected. Use at least one of: --uDOC --sdDOC --medDOC --CovBases --RCN")

    df = pd.read_csv(args.input, sep="\t", comment="#")

    output_cols = {}
    output_cols["locus"] = df["locus"]

    sample_columns = [c for c in df.columns if c != "locus"]

    for sample in sample_columns:

        for field in requested_fields:
            idx = FIELD_MAP[field]
            new_col = f"{sample}_{field}"

            output_cols[new_col] = df[sample].apply(
                lambda x: parse_value(x, idx)
            )

    out_df = pd.DataFrame(output_cols)

    out_df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
