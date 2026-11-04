#!/usr/bin/env python3
"""Region-level circular-shift tests for four prespecified biological processes.

This avoids treating neighboring genes inside the same 160-kb compartment bin
as independent observations. Selection tracks are rotated chromosome-wise,
preserving the number and genomic clustering of observed difference bins.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Source_Data/S2_Table_primary_160kb_compartment_comparison.csv"
GTF = ROOT / "Reference_Data/gencode.v49.annotation.gtf.gz"
REF = ROOT / "Reference_Data"
OUT = ROOT / "Results/prespecified_process_circular_shift.csv"
N_PERM = 20_000
SEED = 20261104


def attrs(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)"', text))


def read_gene_sets() -> dict[str, set[str]]:
    result = {}
    for path in sorted(REF.glob("GO:*.txt")):
        go_id, description = path.stem.split("_", 1)
        result[f"{go_id} {description.replace('_', ' ')}"] = set(path.read_text().split())
    return result


def gene_bins(gene_sets: dict[str, set[str]]) -> dict[str, set[tuple[str, int]]]:
    bins: dict[str, set[tuple[str, int]]] = {name: set() for name in gene_sets}
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            info = attrs(fields[8])
            if info.get("gene_type") != "protein_coding":
                continue
            chrom = fields[0]
            chrom = chrom if chrom.startswith("chr") else f"chr{chrom}"
            if chrom not in {f"chr{i}" for i in range(1, 23)}:
                continue
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            tss = start if strand == "+" else end - 1
            symbol = info.get("gene_name", "")
            bin_start = (tss // 160_000) * 160_000
            for name, symbols in gene_sets.items():
                if symbol in symbols:
                    bins[name].add((chrom, bin_start))
    return bins


def consensus(state: str, positive_fraction: str) -> bool:
    """Return whether a lineage state meets the prespecified 80% consensus rule."""
    if state not in {"A", "B"} or positive_fraction == "":
        return False
    value = float(positive_fraction)
    return value >= 0.8 if state == "A" else value <= 0.2


def tested(row: dict[str, str]) -> bool:
    """Match the analyzable-bin universe used for discovery and gene mapping."""
    return (
        consensus(row["hMSC_state"], row["hMSC_positive_fraction"])
        and consensus(row["AC_state"], row["AC_positive_fraction"])
        and row["AC_vs_hMSC_delta_z"] != ""
    )


def circular_overlap_by_offset(selection: np.ndarray, annotation: np.ndarray) -> np.ndarray:
    # count[k] is the overlap after rotating the selection by k positions.
    raw = np.fft.ifft(np.fft.fft(selection) * np.conj(np.fft.fft(annotation))).real
    return np.rint(raw).astype(int)


def bh_adjust(values: list[float]) -> list[float]:
    n = len(values)
    order = sorted(range(n), key=values.__getitem__)
    adjusted = [1.0] * n
    running = 1.0
    for rank_index in range(n - 1, -1, -1):
        idx = order[rank_index]
        rank = rank_index + 1
        running = min(running, values[idx] * n / rank)
        adjusted[idx] = min(1.0, running)
    return adjusted


def main() -> None:
    gene_sets = read_gene_sets()
    pathway_bins = gene_bins(gene_sets)
    rows_by_chrom: dict[str, list[dict[str, str]]] = defaultdict(list)
    with SOURCE.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["chrom"] in {f"chr{i}" for i in range(1, 23)} and tested(row):
                rows_by_chrom[row["chrom"]].append(row)
    for chrom in rows_by_chrom:
        rows_by_chrom[chrom].sort(key=lambda row: int(row["start"]))

    selections = {
        "all_opposite_state_effect_bins": lambda row: row[
            "AC_vs_hMSC_consensus_effect_difference"
        ].lower() == "true",
        "hMSC_A_to_adipocyte_B": lambda row: row[
            "AC_vs_hMSC_consensus_effect_difference"
        ].lower() == "true" and row["hMSC_state"] == "A" and row["AC_state"] == "B",
        "hMSC_B_to_adipocyte_A": lambda row: row[
            "AC_vs_hMSC_consensus_effect_difference"
        ].lower() == "true" and row["hMSC_state"] == "B" and row["AC_state"] == "A",
    }

    rng = np.random.default_rng(SEED)
    results = []
    for selection_name, predicate in selections.items():
        for pathway_name, bins in pathway_bins.items():
            observed = 0
            null = np.zeros(N_PERM, dtype=int)
            pathway_bin_count = 0
            selected_bin_count = 0
            for chrom, rows in rows_by_chrom.items():
                selection = np.array([predicate(row) for row in rows], dtype=int)
                annotation = np.array(
                    [(chrom, int(row["start"])) in bins for row in rows], dtype=int
                )
                selected_bin_count += int(selection.sum())
                pathway_bin_count += int(annotation.sum())
                counts = circular_overlap_by_offset(selection, annotation)
                observed += int(np.sum(selection * annotation))
                offsets = rng.integers(0, len(rows), size=N_PERM)
                null += counts[offsets]
            p = (1 + int(np.sum(null >= observed))) / (N_PERM + 1)
            null_mean = float(null.mean())
            results.append(
                {
                    "selection": selection_name,
                    "process": pathway_name,
                    "selected_bin_count": selected_bin_count,
                    "process_annotated_bin_count": pathway_bin_count,
                    "observed_overlap_bins": observed,
                    "null_mean_overlap_bins": round(null_mean, 4),
                    "null_2.5_percentile": float(np.quantile(null, 0.025)),
                    "null_97.5_percentile": float(np.quantile(null, 0.975)),
                    "fold_over_null": round(observed / null_mean, 4) if null_mean else "",
                    "empirical_p": p,
                    "permutations": N_PERM,
                    "seed": SEED,
                }
            )
    adjusted = bh_adjust([float(row["empirical_p"]) for row in results])
    for row, value in zip(results, adjusted):
        row["BH_FDR_across_12_tests"] = value
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    for row in results:
        print(
            row["selection"], "|", row["process"], "| observed", row["observed_overlap_bins"],
            "null", row["null_mean_overlap_bins"], "fold", row["fold_over_null"],
            "p", f'{row["empirical_p"]:.4g}', "FDR", f'{row["BH_FDR_across_12_tests"]:.4g}'
        )


if __name__ == "__main__":
    main()
