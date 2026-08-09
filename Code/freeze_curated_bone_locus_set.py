#!/usr/bin/env python3
"""Freeze and audit the 60-autosomal-gene bone-related hypothesis set."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


# Two coordinate-copy errors and two legacy symbols identified during editorial audit.
# Coordinates are GRCh38 midpoint positions from NCBI Gene/Ensembl.
CORRECTIONS = {
    "LRP6": {"official_symbol": "LRP6", "chrom": "chr12", "tss": 12_191_605},
    "MEPE": {"official_symbol": "MEPE", "chrom": "chr4", "tss": 87_834_096},
    "RANKL": {"official_symbol": "TNFSF11"},
    "OPG": {"official_symbol": "TNFRSF11B"},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--loci", required=True, type=Path)
    parser.add_argument("--comparison", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    output = root / "Supplementary/Supplementary_Table_6_curated_bone_loci.csv"
    source = root / "Source_Data/Supplementary_curated_bone_loci_audit.csv"
    output.parent.mkdir(parents=True, exist_ok=True)
    source.parent.mkdir(parents=True, exist_ok=True)

    loci = pd.read_csv(args.loci)
    comparison = pd.read_csv(args.comparison)
    loci = loci.rename(columns={"gene": "gene_symbol_original"})
    loci["official_symbol"] = loci["gene_symbol_original"]
    loci["coordinate_audit"] = "carried_forward_from_frozen_project_set"

    for symbol, correction in CORRECTIONS.items():
        mask = loci["gene_symbol_original"].eq(symbol)
        for column, value in correction.items():
            loci.loc[mask, column] = value
        loci.loc[mask, "coordinate_audit"] = (
            "GRCh38_coordinate_corrected_from_NCBI_Ensembl"
            if "tss" in correction
            else "official_HGNC_symbol_standardized"
        )

    loci["bin_start"] = (loci["tss"].astype(int) // 160_000) * 160_000
    switches = comparison[
        comparison["AC_vs_hMSC_consensus_effect_difference"].astype(bool)
    ][["chrom", "start"]].copy()

    observed = []
    for row in loci.itertuples(index=False):
        local = switches[switches["chrom"].eq(row.chrom)]
        observed.append(bool(((local["start"] - row.bin_start).abs() <= 480_000).any()))
    loci["within_500kb_of_switch"] = observed
    loci["set_definition"] = (
        "curated bone-related hypothesis set assembled from BMD-associated and "
        "monogenic skeletal-disease literature; not a comprehensive GWAS catalog"
    )
    loci["genome_build"] = "GRCh38"
    loci["selection_independence"] = (
        "not independently documented; enrichment is conditional on this frozen set"
    )

    columns = [
        "gene_symbol_original",
        "official_symbol",
        "chrom",
        "tss",
        "bin_start",
        "within_500kb_of_switch",
        "annotation",
        "genome_build",
        "coordinate_audit",
        "set_definition",
        "selection_independence",
    ]
    loci[columns].to_csv(output, index=False)
    loci[columns].to_csv(source, index=False)
    if len(loci) != 60 or int(loci["within_500kb_of_switch"].sum()) != 12:
        raise RuntimeError(
            f"Expected 60 loci and 12 overlaps; observed {len(loci)} and "
            f"{int(loci['within_500kb_of_switch'].sum())}"
        )
    print(output)


if __name__ == "__main__":
    main()
