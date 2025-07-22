#!/usr/bin/env python3
"""
plot_modification_histograms.py

Reads a CSV of triplet modification counts and plots histograms
of total_modifications, stem_modifications, hloop_modifications,
iloop_modifications, bulge_modifications, and mloop_modifications
as subplots in a single PNG. Includes both raw and anchor-sequence-normalized values.
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt

def main():
    parser = argparse.ArgumentParser(
        description="Plot histograms of raw and normalized modification counts from CSV."
    )
    parser.add_argument(
        "input_csv",
        help="Path to the input CSV file (can have commented header lines).",
    )
    parser.add_argument(
        "--output_png",
        default="modification_histograms.png",
        help="Filename for the output PNG (default: %(default)s).",
    )
    args = parser.parse_args()

    # Read CSV, skipping metadata/commented header lines
    df = pd.read_csv(args.input_csv, comment="#")

    # Required columns
    mod_cols = [
        "total_modifications",
        "stem_modifications",
        "hloop_modifications",
        "iloop_modifications",
        "bulge_modifications",
        "mloop_modifications",
    ]

    if "anchor_seq" not in df.columns:
        raise ValueError("Missing required column: anchor_seq")
    missing = [col for col in mod_cols if col not in df.columns]
    if missing:
        raise ValueError(f"Missing expected columns: {missing}")

    # Compute anchor sequence length and normalized modification counts
    df["anchor_len"] = df["anchor_seq"].str.len()
    norm_cols = [f"{col}_norm" for col in mod_cols]
    for col, norm_col in zip(mod_cols, norm_cols):
        df[norm_col] = df[col] / df["anchor_len"]

    # Plot settings
    fig, axes = plt.subplots(nrows=2, ncols=6, figsize=(18, 6))
    axes = axes.flatten()

    # Top row: raw modification counts
    for i, col in enumerate(mod_cols):
        ax = axes[i]
        ax.hist(df[col].dropna(), bins=30, edgecolor="black")
        ax.set_title(f"{col.replace('_', ' ').capitalize()}")
        ax.set_xlabel("Count")
        ax.set_ylabel("Frequency")

    # Bottom row: normalized modification counts
    for i, col in enumerate(norm_cols):
        ax = axes[i + 6]
        ax.hist(df[col].dropna(), bins=30, edgecolor="black")
        ax.set_title(f"{col.replace('_', ' ').capitalize()}")
        ax.set_xlabel("Count / len(anchor_seq)")
        ax.set_ylabel("Frequency")

    plt.tight_layout()
    plt.savefig(args.output_png, dpi=300)
    print(f"Saved histogram PNG to {args.output_png}")

if __name__ == "__main__":
    main()

