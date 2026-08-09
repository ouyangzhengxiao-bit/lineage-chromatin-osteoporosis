#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Figures"
SOURCE = ROOT / "Source_Data"
OUT.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.family": "Arial",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})
def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

COLORS = {"hMSC": "#2F6B9A", "AC": "#D9822B", "OB": "#2E8B57"}

qc = pd.read_csv(SOURCE / "Fig3a_HiC_sample_QC.tsv", sep="\t")
corr = pd.read_csv(SOURCE / "Fig3b_replicate_correlations.csv")
pca = pd.read_csv(SOURCE / "Fig3c_pca_scores.csv")
comp = pd.read_csv(SOURCE / "Fig3d_compartment_comparison.csv")
loci = pd.read_csv(SOURCE / "Fig3e_candidate_locus_comparison.csv")
gwas = pd.read_csv(SOURCE / "Fig3f_gwas_switch_overlap.csv")

sample_order = (
    [f"hMSC_rep{i}" for i in range(1, 8)]
    + [f"AC_rep{i}" for i in range(1, 7)]
    + ["OB_rep1", "OB_rep6"]
)
sample_group = {s: ("hMSC" if s.startswith("hMSC") else "AC" if s.startswith("AC") else "OB")
                for s in sample_order}

mat = pd.DataFrame(np.eye(len(sample_order)), index=sample_order, columns=sample_order)
for row in corr.itertuples(index=False):
    if row.sample_1 in mat.index and row.sample_2 in mat.columns:
        mat.loc[row.sample_1, row.sample_2] = row.spearman_rho
        mat.loc[row.sample_2, row.sample_1] = row.spearman_rho

# Reproduce the chromosome-wise circular-shift null for visualization.
autosomes = [f"chr{i}" for i in range(1, 23)]
valid = comp[comp.chrom.isin(autosomes)].copy()
valid["switch"] = valid.AC_vs_hMSC_consensus_effect_difference.astype(bool)
gwas_auto = gwas[gwas.chrom.isin(autosomes)].copy()
rng = np.random.default_rng(42)
n_perm = 20000
null_hits = np.zeros(n_perm, dtype=int)
for chrom in valid.chrom.drop_duplicates():
    chrom_frame = valid.loc[valid.chrom.eq(chrom)].sort_values("start")
    x = chrom_frame["switch"].to_numpy(dtype=bool)
    if not len(x):
        continue
    expanded = x.copy()
    for offset in (-3, -2, -1, 1, 2, 3):
        expanded |= np.roll(x, offset)
    starts = chrom_frame["start"].to_numpy()
    idx_by_start = {int(v): i for i, v in enumerate(starts)}
    loci_idx = [
        idx_by_start.get(int(v))
        for v in gwas_auto.loc[gwas_auto.chrom.eq(chrom), "bin_start"]
    ]
    loci_idx = [i for i in loci_idx if i is not None]
    if not loci_idx:
        continue
    shifts = rng.integers(0, len(x), size=n_perm)
    shifted_indices = (np.asarray(loci_idx)[None, :] - shifts[:, None]) % len(x)
    null_hits += expanded[shifted_indices].any(axis=1).astype(int)
enrichment_summary = pd.read_csv(SOURCE / "Fig3f_gwas_enrichment.csv").iloc[0]
observed = int(enrichment_summary.observed_hits)
null_mean = float(enrichment_summary.null_mean_hits)
fold = float(enrichment_summary.fold_enrichment)
empirical_p = float(enrichment_summary.empirical_p)

fig = plt.figure(figsize=(11.2, 9.0))
gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.12, 1.0], hspace=0.62, wspace=0.65)
panel_axes = []

# a. Sequencing/contact-map QC
ax = fig.add_subplot(gs[0, :2])
panel_axes.append(ax)
q = qc.copy()
q["group"] = q["sample"].str.split("_").str[0]
q["millions"] = q.contact_sum / 1e6
q["order"] = q["group"].map({"hMSC": 0, "AC": 1, "OB": 2})
q = q.sort_values(["order", "sample"])
bar_colors = [COLORS[g] if inc else "#C7C7C7" for g, inc in zip(q.group, q.included)]
ax.bar(np.arange(len(q)), q.millions, color=bar_colors, width=0.78, edgecolor="none")
ax.axhline(10, color="#777777", lw=0.8, ls="--")
for i, row in enumerate(q.itertuples()):
    if not row.included:
        ax.text(i, max(row.millions, 4), "×", ha="center", va="bottom", color="#8A1C1C",
                fontsize=9, fontweight="bold")
ax.set_ylabel("Valid contacts (millions)")
ax.set_xticks(np.arange(len(q)))
ax.set_xticklabels(q["sample"].str.replace("_rep", " ").tolist(), rotation=60, ha="right")
ax.set_ylim(0, max(q.millions) * 1.12)
despine(ax)

# b. Replicate correlation matrix
ax = fig.add_subplot(gs[0, 2:])
panel_axes.append(ax)
im = ax.imshow(mat.to_numpy(), cmap="magma_r", vmin=0.45, vmax=1.0,
               aspect="equal", interpolation="nearest")
ax.set_xticks(np.arange(len(sample_order)))
ax.set_yticks(np.arange(len(sample_order)))
ax.set_xticklabels([s.replace("_rep", " ") for s in sample_order])
ax.set_yticklabels([s.replace("_rep", " ") for s in sample_order])
cb = fig.colorbar(im, ax=ax, shrink=0.68)
cb.set_label("Spearman ρ")
ax.tick_params(axis="x", rotation=60)
ax.tick_params(axis="y", rotation=0)

# c. PCA
ax = fig.add_subplot(gs[1, :2])
panel_axes.append(ax)
for group, d in pca.groupby("group"):
    ax.scatter(d.PC1, d.PC2, s=48, color=COLORS[group], edgecolor="white",
               linewidth=0.7, label=f"{group} (n={len(d)})", zorder=3)
    if group == "OB":
        for row in d.itertuples():
            ax.text(row.PC1 + 1.6, row.PC2 + 1.3, row.sample.replace("_rep", " "),
                    fontsize=6.2, color=COLORS[group])
ax.axhline(0, color="#D6D6D6", lw=0.6)
ax.axvline(0, color="#D6D6D6", lw=0.6)
ax.set_xlabel("PC1 (45.6%)")
ax.set_ylabel("PC2 (20.3%)")
ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.23))
despine(ax)

# d. Non-inferential hMSC versus adipocyte consensus/effect comparison
ax = fig.add_subplot(gs[1, 2:])
panel_axes.append(ax)
d = comp.dropna(subset=["AC_vs_hMSC_delta_z", "hMSC_positive_fraction", "AC_positive_fraction"]).copy()
d["hMSC_sign_agreement"] = np.maximum(d.hMSC_positive_fraction, 1 - d.hMSC_positive_fraction)
d["AC_sign_agreement"] = np.maximum(d.AC_positive_fraction, 1 - d.AC_positive_fraction)
d["minimum_sign_agreement"] = np.minimum(d.hMSC_sign_agreement, d.AC_sign_agreement)
sw = d.AC_vs_hMSC_consensus_effect_difference.astype(bool)
ax.scatter(d.loc[~sw, "AC_vs_hMSC_delta_z"], d.loc[~sw, "minimum_sign_agreement"],
           s=3, color="#C9CDD2", alpha=0.42, rasterized=True)
ax.scatter(d.loc[sw, "AC_vs_hMSC_delta_z"], d.loc[sw, "minimum_sign_agreement"],
           s=7, color="#A13737", alpha=0.72, rasterized=True)
ax.axvline(-1, color="#777777", ls="--", lw=0.7)
ax.axvline(1, color="#777777", ls="--", lw=0.7)
ax.axhline(0.8, color="#777777", ls="--", lw=0.7)
for gene in ["RSPO3", "VDR", "FAM210A", "CTNNB1"]:
    row = loci[loci.gwas_gene.eq(gene)]
    if len(row) and pd.notna(row.iloc[0].AC_vs_hMSC_delta_z):
        x = row.iloc[0].AC_vs_hMSC_delta_z
        h = max(row.iloc[0].hMSC_positive_fraction, 1 - row.iloc[0].hMSC_positive_fraction)
        a = max(row.iloc[0].AC_positive_fraction, 1 - row.iloc[0].AC_positive_fraction)
        y = min(h, a)
        ax.scatter([x], [y], s=24, color="#111111", zorder=5)
        ax.annotate(gene, (x, y), xytext=(3, 4), textcoords="offset points", fontsize=7)
ax.text(0.98, 0.96, f"{int(sw.sum()):,} consensus differences",
        transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
ax.set_xlabel("Adipocyte − hMSC compartment score (z)")
ax.set_ylabel("Minimum within-lineage sign agreement")
ax.set_ylim(0.45, 1.02)
despine(ax)

# e. Candidate-locus group means across lineages
ax = fig.add_subplot(gs[2, :2])
panel_axes.append(ax)
sel = ["RSPO3", "SATB2", "CTNNB1", "VDR", "FAM210A"]
rows = []
for gene in sel:
    rr = loci[loci.gwas_gene.eq(gene)]
    if rr.empty:
        continue
    key = rr.iloc[0]
    for group, column in (("hMSC", "hMSC_mean_E1"), ("AC", "AC_mean_E1"), ("OB", "OB_mean_E1")):
        rows.append({"gene": gene, "group": group, "value": key[column]})
ld = pd.DataFrame(rows)
pos = {g: i for i, g in enumerate(sel)}
offset = {"hMSC": -0.24, "AC": 0.0, "OB": 0.24}
for group in ["hMSC", "AC", "OB"]:
    dd = ld[ld.group.eq(group)]
    for gene in sel:
        vals = dd.loc[dd.gene.eq(gene), "value"].dropna().to_numpy()
        x0 = pos[gene] + offset[group]
        if len(vals):
            ax.scatter([x0], vals, s=27, color=COLORS[group], alpha=0.86,
                       edgecolor="white", linewidth=0.3)
ax.axhline(0, color="#777777", lw=0.7, ls="--")
ax.set_xticks(range(len(sel)))
ax.set_xticklabels(sel)
ax.set_ylabel("Compartment score (z)")
handles = [mpl.lines.Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS[g],
                           markeredgecolor="none", label=g, markersize=5) for g in ["hMSC", "AC", "OB"]]
ax.legend(handles=handles, frameon=False, ncol=3, loc="upper center",
          bbox_to_anchor=(0.5, -0.18))
despine(ax)

# f. Enrichment null
ax = fig.add_subplot(gs[2, 2:])
panel_axes.append(ax)
bins = np.arange(null_hits.min() - 0.5, max(null_hits.max(), observed) + 1.5, 1)
ax.hist(null_hits, bins=bins, color="#B8C6D1", edgecolor="white", linewidth=0.4)
ax.axvline(observed, color="#A13737", lw=2.0)
ax.text(observed + 0.25, ax.get_ylim()[1] * 0.92,
            f"Exploratory: {observed}/60\n{fold:.2f}-fold (conditional)\nempirical $P = {empirical_p:.2g}$",
        color="#7B1F1F", va="top", fontsize=8)
ax.set_xlabel("Curated loci within ±500 kb of a consensus difference")
ax.set_ylabel("Chromosome-wise shifts")
despine(ax)

for label, ax in zip("ABCDEF", panel_axes):
    ax.text(-0.13, 1.08, label, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top")

fig.savefig(OUT / "Figure3_External_HiC_Validation.png", dpi=300, bbox_inches="tight",
            facecolor="white")
fig.savefig(OUT / "Figure3_External_HiC_Validation.tiff", dpi=300, bbox_inches="tight",
            facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(OUT / "Figure3_External_HiC_Validation.pdf", bbox_inches="tight",
            facecolor="white")
plt.close(fig)

# Submission TIFFs should be explicit RGB, not RGBA.
tiff_path = OUT / "Figure3_External_HiC_Validation.tiff"
with Image.open(tiff_path) as im:
    rgb = Image.new("RGB", im.size, "white")
    if im.mode == "RGBA":
        rgb.paste(im, mask=im.getchannel("A"))
    else:
        rgb.paste(im.convert("RGB"))
    rgb.save(tiff_path, dpi=(300, 300), compression="tiff_lzw")

print(OUT / "Figure3_External_HiC_Validation.png")
