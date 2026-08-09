#!/usr/bin/env python3
"""Generate replicate-aggregated Hi-C maps and compartment tracks at two loci."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import cooler
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd
from PIL import Image


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--hic-root",
        type=Path,
        default=Path(os.environ["HIC_MCOOL_ROOT"]) if "HIC_MCOOL_ROOT" in os.environ else None,
        help="Directory containing SAMPLE/SAMPLE.mcool. May also be HIC_MCOOL_ROOT.",
    )
    parser.add_argument(
        "--e1-table",
        type=Path,
        default=Path(os.environ["HIC_E1_GROUP_TABLE"])
        if "HIC_E1_GROUP_TABLE" in os.environ
        else None,
        help="group_compartment_comparison.csv. May also be HIC_E1_GROUP_TABLE.",
    )
    args = parser.parse_args()
    if args.hic_root is None or args.e1_table is None:
        parser.error("--hic-root and --e1-table (or corresponding environment variables) are required")
    if not args.hic_root.exists():
        parser.error(f"Hi-C root does not exist: {args.hic_root}")
    if not args.e1_table.exists():
        parser.error(f"E1 table does not exist: {args.e1_table}")
    return args


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Figures"
SOURCE = ROOT / "Source_Data"
OUT.mkdir(parents=True, exist_ok=True)
SOURCE.mkdir(parents=True, exist_ok=True)

SAMPLES = {
    "hMSC": [f"hMSC_rep{i}" for i in range(1, 8)],
    "Adipocyte": [f"AC_rep{i}" for i in range(1, 7)],
    "Osteoblast": ["OB_rep1", "OB_rep6"],
}
COLORS = {"hMSC": "#2F6B9A", "Adipocyte": "#D9822B", "Osteoblast": "#2E8B57"}
E1_COLS = {"hMSC": "hMSC_mean_E1", "Adipocyte": "AC_mean_E1", "Osteoblast": "OB_mean_E1"}
LOCI = [
    {"gene": "RSPO3", "chrom": "chr6", "start": 126_200_000, "end": 128_200_000, "tss": 127_159_186},
    {"gene": "VDR", "chrom": "chr12", "start": 46_800_000, "end": 48_800_000, "tss": 47_869_573},
]
CMAP = LinearSegmentedColormap.from_list(
    "hic", ["#FFFFFF", "#F7D6C7", "#E78B6A", "#B33A3A", "#5E1A1A"]
)

plt.rcParams.update(
    {
        "font.family": "Arial",
        "font.size": 8,
        "axes.linewidth": 0.6,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def fetch_balanced(path: Path, locus: dict) -> np.ndarray:
    """Fetch a balanced matrix; never silently substitute unbalanced contacts."""
    clr = cooler.Cooler(f"{path}::/resolutions/20000")
    region = (locus["chrom"], locus["start"], locus["end"])
    try:
        matrix = clr.matrix(balance=True, sparse=False).fetch(region)
    except Exception as error:
        raise RuntimeError(f"Balanced 20-kb matrix unavailable for {path}") from error
    matrix = np.asarray(matrix, dtype=float)
    if not np.isfinite(matrix).any():
        raise RuntimeError(f"Balanced matrix contains no finite contacts: {path}")
    return np.log1p(matrix)


def panel_label(ax, label, x=-0.07, y=1.03):
    ax.text(x, y, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")


def main() -> None:
    args = parse_args()
    e1 = pd.read_csv(args.e1_table)
    sample_paths = {
        group: [args.hic_root / sample / f"{sample}.mcool" for sample in samples]
        for group, samples in SAMPLES.items()
    }
    missing = [str(path) for paths in sample_paths.values() for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing retained mcool files:\n" + "\n".join(missing))

    fig = plt.figure(figsize=(11.5, 7.4))
    outer = fig.add_gridspec(2, 1, hspace=0.30)
    panel_letters = iter("ABCDEFGHIJ")
    matrix_source_rows = []
    group_source_rows = []
    e1_source_rows = []
    balance_rows = []

    for row, locus in enumerate(LOCI):
        sub = outer[row].subgridspec(
            2, 4, height_ratios=[4.0, 1.2], hspace=0.05, wspace=0.12
        )
        replicate_matrices = {}
        group_medians = {}
        starts = np.arange(locus["start"], locus["end"], 20_000)

        for group, paths in sample_paths.items():
            replicate_matrices[group] = []
            for path in paths:
                matrix = fetch_balanced(path, locus)
                replicate_matrices[group].append(matrix)
                balance_rows.append(
                    {
                        "locus": locus["gene"],
                        "sample": path.stem,
                        "group": group,
                        "resolution_bp": 20_000,
                        "balanced": True,
                        "finite_fraction": float(np.isfinite(matrix).mean()),
                    }
                )
                row_index, col_index = np.indices(matrix.shape)
                matrix_source_rows.append(
                    pd.DataFrame(
                        {
                            "locus": locus["gene"],
                            "sample": path.stem,
                            "group": group,
                            "bin1_start": starts[row_index.ravel()],
                            "bin2_start": starts[col_index.ravel()],
                            "log1p_balanced_contact": matrix.ravel(),
                        }
                    )
                )
            group_medians[group] = np.nanmedian(np.stack(replicate_matrices[group]), axis=0)

        difference = group_medians["Adipocyte"] - group_medians["hMSC"]
        positive = np.concatenate(
            [matrix[np.isfinite(matrix) & (matrix > 0)] for matrix in group_medians.values()]
        )
        vmax = np.quantile(positive, 0.995) if len(positive) else 1.0
        diff_limit = np.nanquantile(np.abs(difference), 0.995)
        diff_limit = max(float(diff_limit), 1e-6)

        for col, group in enumerate(("hMSC", "Adipocyte", "Osteoblast")):
            matrix = group_medians[group]
            row_index, col_index = np.indices(matrix.shape)
            group_source_rows.append(
                pd.DataFrame(
                    {
                        "locus": locus["gene"],
                        "matrix": f"{group}_median",
                        "bin1_start": starts[row_index.ravel()],
                        "bin2_start": starts[col_index.ravel()],
                        "value": matrix.ravel(),
                    }
                )
            )
            ax = fig.add_subplot(sub[0, col])
            image = ax.imshow(
                matrix,
                origin="lower",
                cmap=CMAP,
                vmin=0,
                vmax=vmax,
                interpolation="nearest",
                rasterized=True,
                extent=[
                    locus["start"] / 1e6,
                    locus["end"] / 1e6,
                    locus["start"] / 1e6,
                    locus["end"] / 1e6,
                ],
            )
            tss = locus["tss"] / 1e6
            ax.axvline(tss, color="#23303F", lw=0.75, ls="--")
            ax.axhline(tss, color="#23303F", lw=0.75, ls="--")
            ax.text(
                0.03,
                0.96,
                f"{group} median (n={len(SAMPLES[group])})",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.2,
                fontweight="bold",
                bbox=dict(facecolor="white", alpha=0.82, edgecolor="none", pad=1.5),
            )
            panel_label(ax, next(panel_letters))
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_linewidth(0.5)
            if col == 2:
                colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.025)
                colorbar.set_label("log(1 + balanced contact)", fontsize=6.7)
                colorbar.ax.tick_params(labelsize=6)

        row_index, col_index = np.indices(difference.shape)
        group_source_rows.append(
            pd.DataFrame(
                {
                    "locus": locus["gene"],
                    "matrix": "Adipocyte_minus_hMSC",
                    "bin1_start": starts[row_index.ravel()],
                    "bin2_start": starts[col_index.ravel()],
                    "value": difference.ravel(),
                }
            )
        )
        ax = fig.add_subplot(sub[0, 3])
        diff_image = ax.imshow(
            difference,
            origin="lower",
            cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-diff_limit, vcenter=0, vmax=diff_limit),
            interpolation="nearest",
            rasterized=True,
            extent=[
                locus["start"] / 1e6,
                locus["end"] / 1e6,
                locus["start"] / 1e6,
                locus["end"] / 1e6,
            ],
        )
        tss = locus["tss"] / 1e6
        ax.axvline(tss, color="#23303F", lw=0.75, ls="--")
        ax.axhline(tss, color="#23303F", lw=0.75, ls="--")
        ax.text(
            0.03,
            0.96,
            "Adipocyte − hMSC",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.2,
            fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.82, edgecolor="none", pad=1.5),
        )
        panel_label(ax, next(panel_letters))
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.5)
        colorbar = fig.colorbar(diff_image, ax=ax, fraction=0.046, pad=0.025)
        colorbar.set_label("Δ log contact", fontsize=6.7)
        colorbar.ax.tick_params(labelsize=6)

        ax = fig.add_subplot(sub[1, :])
        local = e1[
            (e1["chrom"] == locus["chrom"])
            & (e1["end"] >= locus["start"])
            & (e1["start"] <= locus["end"])
        ].copy()
        x = (local["start"] + local["end"]) / 2 / 1e6
        source_block = local[["chrom", "start", "end"]].copy()
        source_block.insert(0, "locus", locus["gene"])
        for group in COLORS:
            source_block[group] = local[E1_COLS[group]].astype(float).to_numpy()
        e1_source_rows.append(source_block)
        for group, color in COLORS.items():
            y = local[E1_COLS[group]].astype(float)
            ax.plot(x, y, color=color, lw=1.4, label=group)
            ax.fill_between(x, 0, y, color=color, alpha=0.08)
        ax.axhline(0, color="#6B7280", lw=0.6)
        ax.axvline(tss, color="#23303F", lw=0.8, ls="--")
        ax.text(
            tss,
            ax.get_ylim()[1] * 0.88,
            locus["gene"],
            ha="center",
            va="top",
            fontweight="bold",
            fontsize=8,
        )
        panel_label(ax, next(panel_letters), x=-0.012, y=1.06)
        ax.set_xlim(locus["start"] / 1e6, locus["end"] / 1e6)
        ax.set_ylabel("Compartment E1 (z)")
        ax.set_xlabel(f"{locus['chrom']} position (Mb)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        if row == 0:
            ax.legend(
                loc="upper right",
                bbox_to_anchor=(0.995, 0.99),
                frameon=True,
                facecolor="white",
                framealpha=0.88,
                edgecolor="none",
                fontsize=7.2,
                handlelength=1.5,
            )

    fig.subplots_adjust(left=0.065, right=0.97, top=0.97, bottom=0.07)
    pd.concat(matrix_source_rows, ignore_index=True).to_csv(
        SOURCE / "Fig2_replicate_contact_matrices_20kb.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.concat(group_source_rows, ignore_index=True).to_csv(
        SOURCE / "Fig2_group_median_and_difference_matrices_20kb.csv.gz",
        index=False,
        compression="gzip",
    )
    pd.concat(e1_source_rows, ignore_index=True).to_csv(
        SOURCE / "Fig2_compartment_E1_tracks.csv", index=False
    )
    pd.DataFrame(balance_rows).to_csv(SOURCE / "Fig2_balance_audit.csv", index=False)

    png = OUT / "Figure2_HiC_Contact_Maps.png"
    tiff = OUT / "Figure2_HiC_Contact_Maps.tiff"
    pdf = OUT / "Figure2_HiC_Contact_Maps.pdf"
    fig.savefig(png, dpi=300, facecolor="white")
    fig.savefig(tiff, dpi=300, facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(pdf, facecolor="white")
    plt.close(fig)

    with Image.open(tiff) as image:
        image.convert("RGB").save(tiff, dpi=(300, 300), compression="tiff_lzw")
    print(png)


if __name__ == "__main__":
    main()
