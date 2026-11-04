#!/usr/bin/env python3
"""Create figures for the transparent reproducibility-focused manuscript."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT / "Source_Data"
RESULTS = ROOT / "Results"
FIG = ROOT / "Figures"
FIG.mkdir(parents=True, exist_ok=True)

BLUE = "#2B6CB0"
ORANGE = "#DD6B20"
GREEN = "#2F855A"
RED = "#C53030"
GRAY = "#718096"
LIGHT = "#EDF2F7"


def panel(ax, letter: str) -> None:
    ax.text(-0.11, 1.08, letter, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def save(fig, name: str) -> None:
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{name}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(FIG / f"{name}.tiff", dpi=300, bbox_inches="tight", facecolor="white", pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)


def workflow() -> None:
    fig, ax = plt.subplots(figsize=(11, 5.2))
    ax.set_xlim(0, 11); ax.set_ylim(0, 5.2); ax.axis("off")
    boxes = [
        (0.3, 3.4, 2.3, 1.1, "Independent in situ Hi-C\nGSE151319 / SRP260624\n19 libraries; 15 retained", BLUE),
        (3.1, 3.4, 2.2, 1.1, "GC-oriented E1 at 160 kb\nwithin-lineage consensus\nand effect-size filtering", GREEN),
        (5.8, 3.4, 2.1, 1.1, "12,608 tested bins\n777 opposite-state\ndifference bins", ORANGE),
        (8.4, 3.4, 2.2, 1.1, "Genome-wide evaluation\nwithout candidate-gene\npreselection", RED),
        (1.4, 1.1, 2.6, 1.1, "Cross-dataset comparison\nGSE160394, 1-Mb tracks\nhg19→GRCh38", BLUE),
        (4.3, 1.1, 2.5, 1.1, "Complete eBMD atlas\n1,106 signals at 518 loci\nhg19→GRCh38", GREEN),
        (7.1, 1.1, 2.7, 1.1, "Prespecified ontology sets\nskeletal development, ossification,\nfat and mesenchymal differentiation", ORANGE),
    ]
    for x, y, w, h, text, color in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.04,rounding_size=0.08",facecolor="white",edgecolor=color,linewidth=2))
        ax.text(x+w/2,y+h/2,text,ha="center",va="center",fontsize=9,color="#1A202C")
    for x1,y1,x2,y2 in [(2.6,3.95,3.1,3.95),(5.3,3.95,5.8,3.95),(7.9,3.95,8.4,3.95),(9.5,3.4,2.7,2.2),(9.5,3.4,5.55,2.2),(9.5,3.4,8.45,2.2)]:
        ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle="-|>",mutation_scale=12,color=GRAY,linewidth=1.5,connectionstyle="arc3,rad=0"))
    ax.text(5.5,4.95,"Analysis framework and evidence hierarchy",ha="center",fontsize=15,fontweight="bold")
    ax.text(5.5,0.35,"All enrichment and overlap tests used chromosome-wise circular shifts across the analyzable-bin universe (20,000 permutations).",ha="center",fontsize=9,color="#4A5568")
    save(fig, "Figure1_Analysis_Framework")


def technical_validation() -> None:
    qc = pd.read_csv(REPO / "S1_Table_HiC_sample_QC.tsv", sep="\t")
    pairs = pd.read_csv(REPO / "Fig3b_replicate_correlations.csv")
    pca = pd.read_csv(REPO / "Fig3c_pca_scores.csv")
    bins = pd.read_csv(REPO / "S2_Table_primary_160kb_compartment_comparison.csv")
    colors = {"hMSC": BLUE, "AC": ORANGE, "OB": GREEN}
    fig, axes = plt.subplots(2,2,figsize=(10.5,7.4))
    ax=axes[0,0]; panel(ax,"A")
    q=qc.sort_values(["group","sample"]); c=[colors[g] if inc else "#CBD5E0" for g,inc in zip(q.group,q.included)]
    ax.bar(np.arange(len(q)),q.contact_sum/1e6,color=c,edgecolor="white")
    ax.axhline(10,color=RED,ls="--",lw=1); ax.set_ylabel("Valid contacts (millions)"); ax.set_xticks(np.arange(len(q))); ax.set_xticklabels(q["sample"],rotation=70,ha="right",fontsize=7); ax.set_title("Library quality control")
    ax=axes[0,1]; panel(ax,"B")
    same=pairs[pairs.group_1==pairs.group_2].copy(); order=["hMSC","AC","OB"]
    vals=[same.loc[same.group_1==g,"spearman_rho"] for g in order]
    bp=ax.boxplot(vals,tick_labels=["hMSC","Adipocyte","Osteoblast"],patch_artist=True,showfliers=False)
    for patch,g in zip(bp["boxes"],order): patch.set_facecolor(colors[g]); patch.set_alpha(.75)
    ax.set_ylim(0.65,1.01); ax.set_ylabel("Within-lineage Spearman ρ"); ax.set_title("Technical-library concordance")
    ax=axes[1,0]; panel(ax,"C")
    for g,f in pca.groupby("group"):
        ax.scatter(f.PC1,f.PC2,s=45,color=colors[g],label={"AC":"Adipocyte","hMSC":"hMSC","OB":"Osteoblast"}[g],edgecolor="white",linewidth=.5)
    ax.set_xlabel("PC1 (45.6%)"); ax.set_ylabel("PC2 (20.3%)"); ax.legend(frameon=False,fontsize=8); ax.set_title("Lineage separation")
    ax=axes[1,1]; panel(ax,"D")
    tested=bins[bins.AC_vs_hMSC_delta_z.notna()].copy(); sel=tested.AC_vs_hMSC_consensus_effect_difference.astype(bool)
    minagree=np.minimum(np.where(tested.hMSC_state=="A",tested.hMSC_positive_fraction,1-tested.hMSC_positive_fraction),np.where(tested.AC_state=="A",tested.AC_positive_fraction,1-tested.AC_positive_fraction))
    ax.scatter(tested.loc[~sel,"AC_vs_hMSC_delta_z"],minagree[~sel],s=4,color="#CBD5E0",alpha=.5,rasterized=True)
    ax.scatter(tested.loc[sel,"AC_vs_hMSC_delta_z"],minagree[sel],s=8,color=RED,alpha=.8,rasterized=True,label="777 selected bins")
    ax.axvline(-1,color=GRAY,ls="--",lw=1); ax.axvline(1,color=GRAY,ls="--",lw=1); ax.axhline(.8,color=GRAY,ls="--",lw=1)
    ax.set_xlabel("Adipocyte − hMSC compartment Δz"); ax.set_ylabel("Minimum sign agreement"); ax.set_ylim(0,1.03); ax.legend(frameon=False,fontsize=8); ax.set_title("Consensus-and-effect screen")
    fig.tight_layout(); save(fig,"Figure2_Independent_HiC_Reanalysis")


def cross_dataset() -> None:
    comp=pd.read_csv(RESULTS/"cross_dataset_1mb_compartment_comparison.csv")
    null=pd.read_csv(RESULTS/"cross_dataset_overlap_circular_shift_null.csv")["overlap_count"]
    fig,axes=plt.subplots(1,3,figsize=(11,3.5))
    ax=axes[0]; panel(ax,"A")
    ax.hexbin(comp.primary_delta_z_oriented,comp.external_delta_z_mean,gridsize=38,cmap="Blues",mincnt=1)
    ax.axhline(0,color=GRAY,lw=.8); ax.axvline(0,color=GRAY,lw=.8)
    ax.set_xlabel("GSE160394 adipocyte − hMSC Δz"); ax.set_ylabel("GSE151319 adipocyte − hMSC Δz"); ax.set_title("Genome-wide change concordance\nSpearman ρ = −0.070")
    ax=axes[1]; panel(ax,"B")
    vals=[361,107,55]; labels=["Primary switches\nwith coverage","Also switch in\nexternal data","Same switch\ndirection"]
    ax.bar(np.arange(3),vals,color=[BLUE,ORANGE,GREEN]); ax.set_xticks(np.arange(3)); ax.set_xticklabels(labels,fontsize=8); ax.set_ylabel("1-Mb regions"); ax.set_title("Cross-dataset switch replication")
    for i,v in enumerate(vals): ax.text(i,v+8,str(v),ha="center",fontsize=9)
    ax=axes[2]; panel(ax,"C")
    ax.hist(null,bins=np.arange(null.min()-.5,null.max()+1.5,2),color="#CBD5E0",edgecolor="white")
    ax.axvline(111,color=RED,lw=2,label="Observed = 111"); ax.axvline(null.mean(),color=BLUE,ls="--",lw=1.5,label=f"Null mean = {null.mean():.1f}")
    ax.set_xlabel("Overlapping external selected bins"); ax.set_ylabel("Permutations"); ax.set_title("Regional overlap test\nEmpirical P = 0.178"); ax.legend(frameon=False,fontsize=8)
    fig.tight_layout(); save(fig,"Figure3_Cross_Dataset_Reproducibility")


def genetics_and_processes() -> None:
    g=pd.read_csv(RESULTS/"morris_locus_overlap_circular_shift.csv")
    p=pd.read_csv(RESULTS/"prespecified_process_circular_shift.csv")
    p=p[p.selection=="all_opposite_state_effect_bins"].copy()
    p["short"]=["Skeletal system\ndevelopment","Ossification","Fat cell\ndifferentiation","Mesenchymal cell\ndifferentiation"]
    fig,axes=plt.subplots(1,2,figsize=(10.3,4.0))
    ax=axes[0]; panel(ax,"A")
    x=np.arange(2); obs=g.observed_selected_GWAS_bins.to_numpy(); mean=g.null_mean.to_numpy(); lo=mean-g["null_2.5_percentile"].to_numpy(); hi=g["null_97.5_percentile"].to_numpy()-mean
    ax.bar(x-.18,obs,.36,label="Observed",color=RED); ax.bar(x+.18,mean,.36,label="Circular-shift mean",color=BLUE,yerr=np.vstack([lo,hi]),capsize=3)
    ax.set_xticks(x); ax.set_xticklabels(["Signal-containing\n160-kb bins","Within 500 kb of\na signal"],fontsize=8); ax.set_ylabel("Selected compartment-difference bins"); ax.set_title("Complete 518-locus eBMD atlas"); ax.legend(frameon=False,fontsize=8)
    for i,(o,m) in enumerate(zip(obs,mean)): ax.text(i-.18,o+4,str(int(o)),ha="center",fontsize=8); ax.text(i+.18,m+4,f"{m:.1f}",ha="center",fontsize=8)
    ax=axes[1]; panel(ax,"B")
    x=np.arange(len(p)); obs=p.observed_overlap_bins.to_numpy(); mean=p.null_mean_overlap_bins.to_numpy(); lo=mean-p["null_2.5_percentile"].to_numpy(); hi=p["null_97.5_percentile"].to_numpy()-mean
    ax.bar(x-.18,obs,.36,label="Observed",color=ORANGE); ax.bar(x+.18,mean,.36,label="Circular-shift mean",color=GREEN,yerr=np.vstack([lo,hi]),capsize=3)
    ax.set_xticks(x); ax.set_xticklabels(p.short,fontsize=8); ax.set_ylabel("Selected bins overlapping ontology genes"); ax.set_title("Prespecified biological processes\nAll FDR ≥ 0.996"); ax.legend(frameon=False,fontsize=8)
    fig.tight_layout(); save(fig,"Figure4_Genetic_and_Process_Tests")


if __name__ == "__main__":
    plt.rcParams.update({"font.family":"Arial","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42})
    workflow(); technical_validation(); cross_dataset(); genetics_and_processes()
