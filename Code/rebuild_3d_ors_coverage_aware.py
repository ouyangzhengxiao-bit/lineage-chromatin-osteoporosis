#!/usr/bin/env python3
"""Rebuild 3D-ORS with explicit assay availability and no zero-imputation.

The historical score table encoded unavailable scRNA and proteomic measurements as
zero. This script restores missingness, recomputes scRNA specificity from the raw
log-normalized matrix, and reports both a strict five-layer score and an explicitly
coverage-normalized score for candidates with at least four measured layers.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


LAYERS = {
    "L1": ("L1_3d_score", "L1_measured", 0.30),
    "L2": ("L2_gwas_score", "L2_measured", 0.25),
    "L3": ("L3_expr_score", "L3_measured", 0.20),
    "L4": ("L4_chip_score", "L4_measured", 0.15),
    "L5": ("L5_prot_score", "L5_measured", 0.10),
}
BONE_CELL_TYPES = {
    "Osteoblast",
    "Chondrocyte",
    "Adipocyte",
    "Fibroblast",
    "Pericyte",
}


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--h5ad", required=True, type=Path)
    return parser.parse_args()


def compute_l3(data: pd.DataFrame, h5ad: Path) -> pd.DataFrame:
    obj = ad.read_h5ad(h5ad, backed="r")
    if obj.raw is None:
        obj.file.close()
        raise RuntimeError("Processed h5ad lacks raw log-normalized expression")
    cell_type = obj.obs["cell_type"].astype(str).to_numpy()
    available_bone = sorted(BONE_CELL_TYPES.intersection(set(cell_type)))
    other = sorted(set(cell_type).difference(BONE_CELL_TYPES))
    present = [gene for gene in data["gene"] if gene in obj.raw.var_names]
    present_indices = [obj.raw.var_names.get_loc(gene) for gene in present]
    # Slice all requested genes in one pass. Backed CSR column slices materialize
    # the sparse matrix, so looping by gene would repeatedly scan a ~2 GB matrix.
    matrix = obj.raw.X[:, present_indices]
    rows = []
    for gene in data["gene"]:
        if gene not in present:
            rows.append(
                {"gene": gene, "L3_measured": False, "L3_expr_score": np.nan,
                 "L3_bone_expr": np.nan, "L3_other_expr": np.nan,
                 "L3_specificity_ratio": np.nan}
            )
            continue
        idx = present.index(gene)
        vector = matrix[:, idx]
        bone_means = [float(vector[cell_type == ct].mean()) for ct in available_bone]
        other_means = [float(vector[cell_type == ct].mean()) for ct in other]
        mean_bone = float(np.mean(bone_means)) if bone_means else np.nan
        mean_other = float(np.mean(other_means)) if other_means else np.nan
        specificity = mean_bone / (mean_other + 0.01)
        if specificity > 1:
            score = min(np.log2(specificity + 1) / 3, 1.0)
        else:
            score = max(np.log2(specificity + 0.01) / 3, 0.0)
        rows.append(
            {"gene": gene, "L3_measured": True,
             "L3_expr_score": round(float(score), 4),
             "L3_bone_expr": round(mean_bone, 4),
             "L3_other_expr": round(mean_other, 4),
             "L3_specificity_ratio": round(float(specificity), 4)}
        )
    obj.file.close()
    return pd.DataFrame(rows)


def main() -> None:
    a = args()
    root = Path(__file__).resolve().parents[1]
    table = pd.read_csv(a.scores)

    # L1, L2 and L4 were evaluated across the frozen 91-gene candidate universe.
    # L3 and L5 are assay-dependent and therefore retain explicit missing values.
    old_l3 = [c for c in table.columns if c.startswith("L3_")]
    table = table.drop(columns=old_l3, errors="ignore")
    table = table.merge(compute_l3(table, a.h5ad), on="gene", how="left")
    table["L1_measured"] = table["L1_3d_score"].notna()
    table["L2_measured"] = table["L2_gwas_score"].notna()
    table["L4_measured"] = table["L4_chip_score"].notna()
    table["L5_measured"] = table["L5_in_proteomics"].fillna(False).astype(bool)
    for column in ["L5_prot_score", "L5_BMD_logFC", "L5_BMD_P"]:
        table.loc[~table["L5_measured"], column] = np.nan

    for layer, (score, measured, weight) in LAYERS.items():
        table[f"{layer}_contrib"] = np.where(
            table[measured], table[score] * weight, np.nan
        )
    measured_cols = [v[1] for v in LAYERS.values()]
    table["evidence_layers_measured"] = table[measured_cols].sum(axis=1).astype(int)
    table["measured_layers"] = table.apply(
        lambda row: ";".join(layer for layer, (_, flag, _) in LAYERS.items() if row[flag]),
        axis=1,
    )
    weighted_sum = sum(
        table[score].fillna(0) * weight for score, _, weight in LAYERS.values()
    )
    measured_weight = sum(
        table[flag].astype(float) * weight for _, flag, weight in LAYERS.values()
    )
    table["3D_ORS_available"] = np.where(
        table["evidence_layers_measured"] >= 4, weighted_sum / measured_weight, np.nan
    )
    table["3D_ORS_strict_5layer"] = np.where(
        table["evidence_layers_measured"] == 5, weighted_sum, np.nan
    )
    table["candidate_source"] = "predefined 3D-compartment candidate universe"
    table["candidate_universe_note"] = (
        "91 genes carried forward from the compartment-screen pipeline; not an independent validation set"
    )
    table["within_coverage_group_rank"] = pd.Series(pd.NA, index=table.index, dtype="Int64")
    for layers in (4, 5):
        group_mask = (
            table["evidence_layers_measured"].eq(layers)
            & table["3D_ORS_available"].notna()
        )
        table.loc[group_mask, "within_coverage_group_rank"] = table.loc[
            group_mask, "3D_ORS_available"
        ].rank(method="min", ascending=False).astype("Int64")
    table["ranking_scope"] = table["evidence_layers_measured"].map(
        {3: "not ranked (<4 measured layers)",
         4: "ranked within four-layer candidates",
         5: "ranked within five-layer candidates"}
    )
    strict_mask = table["3D_ORS_strict_5layer"].notna()
    table["strict_rank"] = pd.Series(pd.NA, index=table.index, dtype="Int64")
    table.loc[strict_mask, "strict_rank"] = table.loc[
        strict_mask, "3D_ORS_strict_5layer"
    ].rank(method="min", ascending=False).astype("Int64")

    # Retain the historical column only as an audit trail; it is not used downstream.
    if "3D_ORS" in table:
        table = table.rename(columns={"3D_ORS": "historical_3D_ORS_zero_imputed"})
    if "ors_quartile" in table:
        table = table.rename(columns={"ors_quartile": "historical_ors_quartile"})

    ordered = [
        "gene", "candidate_source", "candidate_universe_note",
        "evidence_layers_measured", "measured_layers", "3D_ORS_available",
        "within_coverage_group_rank", "ranking_scope",
        "3D_ORS_strict_5layer", "strict_rank",
    ]
    remaining = [c for c in table.columns if c not in ordered]
    table = table[ordered + remaining].sort_values(
        ["3D_ORS_available", "gene"], ascending=[False, True], na_position="last"
    )
    table.to_csv(root / "Supplementary/Supplementary_Table_2_3D_ORS_scores.csv", index=False)
    table.to_csv(root / "Source_Data/Fig5_3D_ORS_scores.csv", index=False)

    counts = table["evidence_layers_measured"].value_counts().sort_index()
    counts.rename_axis("evidence_layers_measured").rename("n_candidates").reset_index().to_csv(
        root / "Source_Data/Fig5_measurement_coverage.csv", index=False
    )
    print(f"Candidates: {len(table)}")
    print(f"Coverage-aware (>=4 layers): {table['3D_ORS_available'].notna().sum()}")
    print(f"Strict five-layer: {table['3D_ORS_strict_5layer'].notna().sum()}")
    print(table[["gene", "evidence_layers_measured", "3D_ORS_available", "3D_ORS_strict_5layer"]].head(15).to_string(index=False))


if __name__ == "__main__":
    main()
