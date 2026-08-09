#!/usr/bin/env python3
"""Summarize GC-oriented 160 kb Hi-C technical-library comparisons."""

from __future__ import annotations

import glob
import itertools
import os
from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
WORK_ROOT = Path(os.environ.get("EXTERNAL_HIC_WORK", "work"))
EIG_DIR = WORK_ROOT / "external_validation_160kb/eigs"
OUT_DIR = WORK_ROOT / "external_validation_160kb/analysis"
GWAS_PATH = PACKAGE_ROOT / "Source_Data/Fig3e_candidate_locus_comparison.csv"
CURATED_LOCI_PATH = PACKAGE_ROOT / "Supplementary/Supplementary_Table_6_curated_bone_loci.csv"


def consensus_state(values: np.ndarray, threshold: float = 0.8) -> tuple[np.ndarray, np.ndarray]:
    finite = np.isfinite(values)
    n = finite.sum(axis=0)
    pos_fraction = np.divide(
        (values > 0).sum(axis=0),
        n,
        out=np.full(values.shape[1], np.nan),
        where=n > 0,
    )
    state = np.where(pos_fraction >= threshold, "A", np.where(pos_fraction <= 1 - threshold, "B", "ambiguous"))
    state[n == 0] = "missing"
    return state, pos_fraction


def load_curated_gwas_loci() -> dict:
    table = pd.read_csv(CURATED_LOCI_PATH)
    return {
        row.official_symbol: (row.chrom, int(row.tss), int(row.tss), row.annotation)
        for row in table.itertuples(index=False)
    }


def circular_shift_enrichment(
    comparison: pd.DataFrame,
    gwas_loci: dict,
    switch_column: str,
    permutations: int = 20_000,
    seed: int = 42,
) -> tuple[pd.DataFrame, dict]:
    loci_rows = []
    by_chromosome = {}
    for chromosome, frame in comparison.groupby("chrom", sort=False):
        ordered = frame.sort_values("start").reset_index(drop=True)
        switch = ordered[switch_column].astype(bool).to_numpy()
        expanded = switch.copy()
        for offset in (-3, -2, -1, 1, 2, 3):
            expanded |= np.roll(switch, offset)
        by_chromosome[chromosome] = {
            "n": len(ordered),
            "start_to_index": {int(start): index for index, start in enumerate(ordered["start"])},
            "expanded": expanded,
        }

    gene_indices = {}
    for gene, (chromosome, start, end, note) in gwas_loci.items():
        if chromosome not in by_chromosome:
            continue
        tss = (int(start) + int(end)) // 2
        bin_start = (tss // 160_000) * 160_000
        index = by_chromosome[chromosome]["start_to_index"].get(bin_start)
        if index is None:
            continue
        hit = bool(by_chromosome[chromosome]["expanded"][index])
        loci_rows.append(
            {
                "gene": gene,
                "chrom": chromosome,
                "tss": tss,
                "bin_start": bin_start,
                "within_500kb_of_switch": hit,
                "annotation": note,
            }
        )
        gene_indices.setdefault(chromosome, []).append(index)

    observed = int(sum(row["within_500kb_of_switch"] for row in loci_rows))
    rng = np.random.default_rng(seed)
    null = np.zeros(permutations, dtype=int)
    for chromosome, indices in gene_indices.items():
        indices = np.asarray(indices, dtype=int)
        info = by_chromosome[chromosome]
        shifts = rng.integers(0, info["n"], size=permutations)
        shifted_indices = (indices[None, :] - shifts[:, None]) % info["n"]
        null += info["expanded"][shifted_indices].any(axis=1).astype(int)

    empirical_p = float((1 + np.sum(null >= observed)) / (permutations + 1))
    expected = float(null.mean())
    summary = {
        "n_loci": len(loci_rows),
        "observed_hits": observed,
        "null_mean_hits": expected,
        "fold_enrichment": float(observed / expected) if expected else np.nan,
        "empirical_p": empirical_p,
        "permutations": permutations,
        "window": "±500 kb (three 160 kb bins on each side)",
        "null_model": "independent chromosome-wise circular shifts of the switch track",
    }
    return pd.DataFrame(loci_rows), summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    files = sorted(glob.glob(str(EIG_DIR / "*.cis.vecs.tsv")))
    series = {}
    coords = None
    for file in files:
        sample = os.path.basename(file).split(".")[0]
        frame = pd.read_csv(file, sep="\t")
        frame = frame[frame["chrom"].str.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])")].copy()
        index = pd.MultiIndex.from_frame(frame[["chrom", "start", "end"]])
        series[sample] = pd.Series(frame["E1"].to_numpy(), index=index, name=sample)
        if coords is None:
            coords = index

    raw = pd.concat(series.values(), axis=1)
    raw.columns = list(series)

    # Standardize within chromosome and sample while retaining GC-oriented sign.
    z = raw.copy()
    chrom = raw.index.get_level_values(0)
    for sample in z.columns:
        for chromosome in chrom.unique():
            mask = chrom == chromosome
            values = z.loc[mask, sample]
            mean = values.mean()
            sd = values.std()
            z.loc[mask, sample] = (values - mean) / sd if sd and np.isfinite(sd) else np.nan

    raw.reset_index().to_csv(OUT_DIR / "oriented_E1_by_sample.csv", index=False)
    z.reset_index().to_csv(OUT_DIR / "standardized_E1_by_sample.csv", index=False)

    pair_rows = []
    for first, second in itertools.combinations(raw.columns, 2):
        shared = raw[[first, second]].dropna()
        rho = shared[first].corr(shared[second], method="spearman")
        pair_rows.append(
            {
                "sample_1": first,
                "sample_2": second,
                "group_1": first.split("_", 1)[0],
                "group_2": second.split("_", 1)[0],
                "n_bins": len(shared),
                "spearman_rho": rho,
            }
        )
    pairs = pd.DataFrame(pair_rows)
    pairs.to_csv(OUT_DIR / "pairwise_replicate_correlations.csv", index=False)

    groups = {
        group: [column for column in raw.columns if column.startswith(f"{group}_")]
        for group in ("hMSC", "AC", "OB")
    }
    result = raw.reset_index()[["chrom", "start", "end"]].copy()
    for group, columns in groups.items():
        array = raw[columns].to_numpy().T
        state, pos_fraction = consensus_state(array)
        result[f"{group}_mean_E1"] = np.nanmean(array, axis=0)
        result[f"{group}_state"] = state
        result[f"{group}_positive_fraction"] = pos_fraction
        result[f"{group}_n"] = np.isfinite(array).sum(axis=0)

    for group in ("AC", "OB"):
        control = z[groups["hMSC"]].to_numpy().T
        case = z[groups[group]].to_numpy().T
        result[f"{group}_vs_hMSC_delta_z"] = np.nanmean(case, axis=0) - np.nanmean(control, axis=0)
        opposite = (
            result[f"{group}_state"].isin(["A", "B"])
            & result["hMSC_state"].isin(["A", "B"])
            & (result[f"{group}_state"] != result["hMSC_state"])
        )
        # These are technical libraries, not independent donors. Selection therefore
        # uses only sign consensus and a prespecified effect-size threshold; no P or q
        # values are calculated from technical-library variation.
        result[f"{group}_vs_hMSC_consensus_effect_difference"] = (
            opposite
            & (result[f"{group}_vs_hMSC_delta_z"].abs() >= 1.0)
            & ((result[f"{group}_n"] == 2) if group == "OB" else True)
        )

    result.to_csv(OUT_DIR / "group_compartment_comparison.csv", index=False)

    external_loci, enrichment = circular_shift_enrichment(
        result,
        load_curated_gwas_loci(),
        "AC_vs_hMSC_consensus_effect_difference",
    )
    external_loci.to_csv(OUT_DIR / "curated_gwas_loci_external_switch_overlap.csv", index=False)
    pd.DataFrame([enrichment]).to_csv(
        OUT_DIR / "curated_gwas_circular_shift_enrichment.csv", index=False
    )

    # PCA across common, variable autosomal bins.
    pca_matrix = z.T
    common = pca_matrix.notna().mean(axis=0) >= 0.9
    variable = pca_matrix.loc[:, common].var(axis=0) > 0.05
    x = pca_matrix.loc[:, common].loc[:, variable].fillna(0).to_numpy().copy()
    x -= x.mean(axis=0, keepdims=True)
    u, singular, _ = np.linalg.svd(x, full_matrices=False)
    scores = u[:, :3] * singular[:3]
    explained = singular**2 / np.sum(singular**2)
    pca = pd.DataFrame(
        {
            "sample": pca_matrix.index,
            "group": [sample.split("_", 1)[0] for sample in pca_matrix.index],
            "PC1": scores[:, 0],
            "PC2": scores[:, 1],
            "PC3": scores[:, 2],
        }
    )
    pca.to_csv(OUT_DIR / "sample_pca_scores.csv", index=False)

    gwas = pd.read_csv(GWAS_PATH).loc[:, :"MSC_expr"].drop_duplicates("gwas_gene")
    gwas["start"] = (gwas["gwas_tss"] // 160_000) * 160_000
    gwas["end"] = gwas["start"] + 160_000
    loci = gwas.merge(
        result,
        left_on=["gwas_chrom", "start", "end"],
        right_on=["chrom", "start", "end"],
        how="left",
    )
    loci.to_csv(OUT_DIR / "candidate_locus_comparison.csv", index=False)

    summary = pairs.copy()
    summary[["group_low", "group_high"]] = np.sort(
        summary[["group_1", "group_2"]].to_numpy(), axis=1
    )
    corr_summary = (
        summary.groupby(["group_low", "group_high"])["spearman_rho"]
        .agg(["count", "mean", "median", "min", "max"])
        .reset_index()
    )
    corr_summary.to_csv(OUT_DIR / "correlation_summary.csv", index=False)

    with open(OUT_DIR / "analysis_summary.txt", "w") as handle:
        handle.write("External GSE151319/SRP260624 Hi-C technical comparison at 160 kb\n")
        handle.write(f"Samples analyzed: {len(raw.columns)} ({', '.join(f'{k}={len(v)}' for k, v in groups.items())})\n")
        handle.write(
            f"PCA explained variance: PC1={explained[0]:.3f}, PC2={explained[1]:.3f}, PC3={explained[2]:.3f}\n"
        )
        handle.write(
            "AC versus hMSC consensus-and-effect differences: "
            f"{int(result['AC_vs_hMSC_consensus_effect_difference'].sum())}\n"
        )
        handle.write(
            "Exploratory OB versus hMSC consensus-and-effect differences: "
            f"{int(result['OB_vs_hMSC_consensus_effect_difference'].sum())}\n"
        )
        handle.write(
            "Chromosome-wise circular-shift enrichment: "
            f"{enrichment['observed_hits']}/{enrichment['n_loci']} loci observed; "
            f"null mean={enrichment['null_mean_hits']:.2f}; "
            f"fold={enrichment['fold_enrichment']:.2f}; "
            f"empirical P={enrichment['empirical_p']:.4g}\n"
        )
        handle.write("\nPairwise correlation summary:\n")
        handle.write(corr_summary.to_string(index=False))
        handle.write("\n\nCandidate loci:\n")
        columns = [
            "gwas_gene",
            "hMSC_state",
            "AC_state",
            "OB_state",
            "AC_vs_hMSC_delta_z",
            "AC_vs_hMSC_consensus_effect_difference",
            "OB_vs_hMSC_delta_z",
            "OB_vs_hMSC_consensus_effect_difference",
        ]
        handle.write(loci[columns].to_string(index=False))
        handle.write("\n")


if __name__ == "__main__":
    main()
