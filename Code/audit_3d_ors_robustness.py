#!/usr/bin/env python3
"""Recalculate leave-one-layer-out ranks for strict five-layer candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


LAYERS = {
    "L1_3D_chromatin": ("L1_3d_score", 0.30),
    "L2_genetic_annotation": ("L2_gwas_score", 0.25),
    "L3_scRNA": ("L3_expr_score", 0.20),
    "L4_ChIP": ("L4_chip_score", 0.15),
    "L5_proteomics": ("L5_prot_score", 0.10),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", required=True, type=Path, help="3d_ors_full_scores.csv")
    return parser.parse_args()


def rank_descending(values: pd.Series) -> pd.Series:
    return values.rank(method="min", ascending=False).astype(int)


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    source_dir = root / "Source_Data"
    supplementary_dir = root / "Supplementary"
    source_dir.mkdir(parents=True, exist_ok=True)
    supplementary_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(args.scores)
    required = [column for column, _ in LAYERS.values()]
    complete = data.dropna(subset=required).copy()
    if len(complete) < 5:
        raise RuntimeError(f"Too few strict five-layer candidates: {len(complete)}")

    weights = {name: weight for name, (_, weight) in LAYERS.items()}
    complete["baseline_score_recomputed"] = sum(
        complete[column] * weight for column, weight in LAYERS.values()
    )
    complete["baseline_rank"] = rank_descending(complete["baseline_score_recomputed"])

    summary_rows = []
    gene_rows = []
    baseline_order = complete.sort_values("baseline_rank")
    baseline_top5 = set(baseline_order.head(5)["gene"])
    baseline_top10 = set(baseline_order.head(10)["gene"])

    for omitted, (omitted_column, _) in LAYERS.items():
        retained = {
            name: (column, weight)
            for name, (column, weight) in LAYERS.items()
            if name != omitted
        }
        denominator = sum(weight for _, weight in retained.values())
        score_column = f"score_without_{omitted}"
        rank_column = f"rank_without_{omitted}"
        complete[score_column] = sum(
            complete[column] * (weight / denominator)
            for column, weight in retained.values()
        )
        complete[rank_column] = rank_descending(complete[score_column])
        rho = float(spearmanr(complete["baseline_rank"], complete[rank_column]).statistic)
        ranked = complete.sort_values(rank_column)
        top5 = set(ranked.head(5)["gene"])
        top10 = set(ranked.head(10)["gene"])
        summary_rows.append(
            {
                "omitted_layer": omitted,
                "n_strict_five_layer_candidates": len(complete),
                "spearman_rank_correlation": rho,
                "top5_overlap_n": len(baseline_top5 & top5),
                "top10_overlap_n": len(baseline_top10 & top10),
                "top5_genes": ";".join(ranked.head(5)["gene"]),
                "top10_genes": ";".join(ranked.head(10)["gene"]),
            }
        )
        for row in complete.itertuples(index=False):
            gene_rows.append(
                {
                    "gene": row.gene,
                    "omitted_layer": omitted,
                    "baseline_rank": getattr(row, "baseline_rank"),
                    "leave_one_layer_out_rank": getattr(row, rank_column),
                    "rank_change": getattr(row, rank_column) - getattr(row, "baseline_rank"),
                    "leave_one_layer_out_score": getattr(row, score_column),
                }
            )

    summary = pd.DataFrame(summary_rows)
    gene_level = pd.DataFrame(gene_rows)
    summary.to_csv(
        supplementary_dir / "Supplementary_Table_7_3D_ORS_leave_one_layer_out.csv",
        index=False,
    )
    gene_level.to_csv(
        source_dir / "Supplementary_3D_ORS_leave_one_layer_out_gene_ranks.csv",
        index=False,
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
