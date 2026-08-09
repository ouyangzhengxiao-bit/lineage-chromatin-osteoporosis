#!/usr/bin/env python3
"""Generate the single-cell main figure from the processed GSE317068 h5ad."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from PIL import Image


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=Path(os.environ["GSE317068_H5AD"])
        if "GSE317068_H5AD" in os.environ
        else None,
        help="Processed GSE317068 h5ad. May also be set with GSE317068_H5AD.",
    )
    args = parser.parse_args()
    if args.h5ad is None:
        parser.error("--h5ad or GSE317068_H5AD is required")
    if not args.h5ad.exists():
        parser.error(f"h5ad does not exist: {args.h5ad}")
    return args


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Figures"
SOURCE = ROOT / "Source_Data"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)

CELL_COLORS = {
    "Adipocyte": "#E69F00",
    "Chondrocyte": "#F0E442",
    "Endothelial": "#009E73",
    "Erythroid": "#D55E00",
    "Fibroblast": "#7B2CBF",
    "Lymphocyte": "#8C564B",
    "Myeloid": "#CC79A7",
    "Osteoblast": "#56B4E9",
    "Pericyte": "#00A6B4",
}
CONDITION_COLORS = {
    "Disease": "#2F6B9A",
    "Normal": "#D9822B",
    "Sorted_ALPL": "#2E8B57",
}
EXPR_CMAP = LinearSegmentedColormap.from_list(
    "expression", ["#E8EDF2", "#8BC4D8", "#2E8B57", "#F0C808", "#B33A3A"]
)
FEATURE_GENES = ["RSPO3", "VDR", "TGFBR1", "SMAD3"]
DOT_GENES = ["RUNX2", "ALPL", "BGLAP", "RSPO3", "VDR", "TGFBR1", "SMAD3", "ZBTB40"]

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8.2,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def style_umap(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.015)


def panel_label(ax, label, x=-0.05, y=1.03):
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
    )


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.h5ad, backed="r")
    requested = list(dict.fromkeys(FEATURE_GENES + DOT_GENES))
    missing = [gene for gene in requested if gene not in adata.var_names]
    if missing:
        adata.file.close()
        raise RuntimeError(f"Required genes are missing from the processed object: {missing}")
    if adata.raw is None or "counts" not in adata.layers:
        adata.file.close()
        raise RuntimeError("The h5ad must contain raw log-normalized values and a counts layer")

    umap = np.asarray(adata.obsm["X_umap"])
    cell_type = adata.obs["cell_type"].astype(str).to_numpy()
    condition = adata.obs["condition"].astype(str).to_numpy()
    indices = [adata.var_names.get_loc(gene) for gene in requested]
    lognorm_matrix = adata.raw.X[:, indices]
    count_matrix = adata.layers["counts"][:, indices]
    if hasattr(lognorm_matrix, "toarray"):
        lognorm_matrix = lognorm_matrix.toarray()
    if hasattr(count_matrix, "toarray"):
        count_matrix = count_matrix.toarray()
    lognorm_matrix = np.asarray(lognorm_matrix, dtype=np.float32)
    count_matrix = np.asarray(count_matrix)
    values = {gene: lognorm_matrix[:, i] for i, gene in enumerate(requested)}
    detected = {gene: count_matrix[:, i] > 0 for i, gene in enumerate(requested)}
    adata.file.close()

    embedding = pd.DataFrame(
        {
            "UMAP1": umap[:, 0],
            "UMAP2": umap[:, 1],
            "cell_type": cell_type,
            "sampling_stratum": condition,
            **values,
        }
    )
    embedding.to_csv(
        SOURCE / "Fig1_single_cell_embedding_expression.csv.gz",
        index=False,
        compression="gzip",
    )

    dot_rows = []
    for label in CELL_COLORS:
        mask = cell_type == label
        for gene in DOT_GENES:
            x = values[gene][mask]
            dot_rows.append(
                {
                    "cell_type": label,
                    "gene": gene,
                    "mean_normalized_expression": float(np.mean(x)),
                    "percent_expressing": float(100 * np.mean(detected[gene][mask])),
                    "n_cells": int(mask.sum()),
                }
            )
    dot = pd.DataFrame(dot_rows)
    dot.to_csv(SOURCE / "Fig1_candidate_gene_dotplot.csv", index=False)
    for gene in DOT_GENES:
        mask = dot["gene"].eq(gene)
        x = dot.loc[mask, "mean_normalized_expression"]
        sd = x.std(ddof=0)
        dot.loc[mask, "scaled_mean_expression"] = (
            (x - x.mean()) / sd if sd > 0 else 0.0
        )

    fig = plt.figure(figsize=(11.4, 7.25))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.12, 0.88], hspace=0.12)
    top = outer[0].subgridspec(1, 3, width_ratios=[1.18, 1.18, 1.0], wspace=0.13)
    bottom = outer[1].subgridspec(1, 4, wspace=0.08)

    # A: broad cell-class UMAP.
    ax = fig.add_subplot(top[0, 0])
    for label, color in CELL_COLORS.items():
        idx = np.where(cell_type == label)[0]
        ax.scatter(
            umap[idx, 0],
            umap[idx, 1],
            s=0.45,
            alpha=0.65,
            color=color,
            linewidths=0,
            rasterized=True,
            label=label,
        )
    style_umap(ax)
    panel_label(ax, "A")
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        frameon=True,
        facecolor="white",
        framealpha=0.88,
        edgecolor="none",
        markerscale=7,
        handletextpad=0.3,
        borderaxespad=0,
        ncol=2,
        columnspacing=0.7,
        fontsize=6.7,
    )

    # B: stored sampling strata.
    ax = fig.add_subplot(top[0, 1])
    for label, color in CONDITION_COLORS.items():
        idx = np.where(condition == label)[0]
        ax.scatter(
            umap[idx, 0],
            umap[idx, 1],
            s=0.45,
            alpha=0.55,
            color=color,
            linewidths=0,
            rasterized=True,
            label=label.replace("_", " "),
        )
    style_umap(ax)
    panel_label(ax, "B")
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.01, 0.99),
        frameon=True,
        facecolor="white",
        framealpha=0.88,
        edgecolor="none",
        markerscale=7,
        handletextpad=0.3,
        borderaxespad=0,
        fontsize=7.1,
    )

    # C: candidate-gene expression by broad cell class.
    ax = fig.add_subplot(top[0, 2])
    cell_order = list(CELL_COLORS)
    for _, row in dot.iterrows():
        x = DOT_GENES.index(row["gene"])
        y = cell_order.index(row["cell_type"])
        ax.scatter(
            x,
            y,
            s=8 + row["percent_expressing"] * 1.15,
            c=row["scaled_mean_expression"],
            cmap="RdYlBu_r",
            vmin=-2,
            vmax=2,
            edgecolor="#4A5560",
            linewidth=0.25,
        )
    ax.set_xticks(range(len(DOT_GENES)))
    ax.set_xticklabels(DOT_GENES, rotation=45, ha="left")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", top=True, bottom=False, labeltop=True, labelbottom=False, pad=2)
    ax.set_yticks(range(len(cell_order)))
    ax.set_yticklabels(cell_order)
    ax.invert_yaxis()
    ax.grid(color="#E5E7EB", linewidth=0.45)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_label(ax, "C", x=-0.15)
    size_handles = [
        ax.scatter([], [], s=8 + pct * 1.15, facecolor="#B7C7D4", edgecolor="#4A5560",
                   linewidth=0.25, label=f"{pct}%")
        for pct in (10, 40, 70)
    ]
    ax.legend(
        handles=size_handles,
        title="Expressing",
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.0, 1.0),
        fontsize=6.5,
        title_fontsize=6.5,
        borderaxespad=0,
    )

    # D–G: prioritized genes on the shared embedding.
    for index, gene in enumerate(FEATURE_GENES):
        ax = fig.add_subplot(bottom[0, index])
        x = values[gene]
        finite = np.isfinite(x)
        vmax = np.nanpercentile(x[finite], 99.5)
        order = np.argsort(x, kind="stable")
        sc = ax.scatter(
            umap[order, 0],
            umap[order, 1],
            c=x[order],
            s=0.34,
            cmap=EXPR_CMAP,
            vmin=0,
            vmax=max(vmax, 1e-6),
            linewidths=0,
            rasterized=True,
        )
        style_umap(ax)
        ax.text(
            0.03,
            0.96,
            gene,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=9,
            fontweight="bold",
            color="#1F2933",
        )
        panel_label(ax, chr(ord("D") + index))
        cb = fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.01)
        cb.set_ticks([])
        cb.outline.set_linewidth(0.4)

    fig.subplots_adjust(left=0.035, right=0.95, top=0.86, bottom=0.06)
    png = OUT / "Figure1_SingleCell_Landscape.png"
    tiff = OUT / "Figure1_SingleCell_Landscape.tiff"
    pdf = OUT / "Figure1_SingleCell_Landscape.pdf"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(tiff, dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    with Image.open(tiff) as image:
        image.convert("RGB").save(tiff, dpi=(300, 300), compression="tiff_lzw")
    # The UMAP point layers are raster data; embed the reviewed 300-dpi rendering
    # in PDF rather than creating an unnecessarily large vector point stream.
    with Image.open(png) as image:
        image.convert("RGB").save(pdf, "PDF", resolution=300)
    print(png)


if __name__ == "__main__":
    main()
