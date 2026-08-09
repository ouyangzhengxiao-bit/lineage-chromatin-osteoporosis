#!/usr/bin/env python3
"""Export compact source-data tables from the stored single-cell objects."""

import argparse
import os
from pathlib import Path

import pandas as pd
import scanpy as sc

parser = argparse.ArgumentParser()
parser.add_argument(
    "--project-root",
    type=Path,
    default=Path(os.environ["PROJECT_3D_OSTEOPOROSIS"])
    if "PROJECT_3D_OSTEOPOROSIS" in os.environ
    else None,
)
args = parser.parse_args()
if args.project_root is None or not args.project_root.exists():
    parser.error("--project-root or PROJECT_3D_OSTEOPOROSIS is required")
PROJECT = args.project_root
OUT = Path(__file__).resolve().parents[1] / "Source_Data"
OUT.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(PROJECT / "results/single_cell_analysis/gse317068_processed.h5ad", backed="r")

for column, filename in (
    ("cell_type", "Fig1a_cell_type_counts.csv"),
    ("condition", "Fig1b_condition_counts.csv"),
    ("donor", "Supplementary_donor_counts.csv"),
    ("site", "Supplementary_site_counts.csv"),
):
    counts = adata.obs[column].value_counts(dropna=False).rename_axis(column).reset_index(name="n_cells")
    counts.to_csv(OUT / filename, index=False)

degs = pd.read_csv(PROJECT / "results/single_cell_analysis/disease_vs_normal_DEGs.csv")
degs.to_csv(OUT / "Fig1c_disease_vs_normal_DEGs.csv", index=False)

switches = pd.read_csv(PROJECT / "results/3d_genome_integration/compartment_switch_genes.csv")
summary = (
    switches.groupby(["comparison", "switch_type"], dropna=False)
    .size()
    .rename("n_gene_assignments")
    .reset_index()
)
summary.to_csv(OUT / "Fig1d_switch_gene_assignments.csv", index=False)

print(f"Source-data tables written to {OUT}")
