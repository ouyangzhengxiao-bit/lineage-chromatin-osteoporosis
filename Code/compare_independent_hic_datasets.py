#!/usr/bin/env python3
"""Compare adipocyte-versus-stromal compartment changes across two Hi-C studies.

GSE160394 provides one processed 1-Mb compartment track per lineage in hg19.
Those intervals are lifted to GRCh38 and compared with the independently
processed 160-kb GSE151319/SRP260624 consensus map. The comparison is explicitly
regional and descriptive because the primary study has one processed track per
lineage and the external libraries are technical, not donor, replicates.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
from pyliftover import LiftOver


PRIMARY = ROOT / "Reference_Data/GSE160394_processed"
EXTERNAL = ROOT / "Source_Data/S2_Table_primary_160kb_compartment_comparison.csv"
CHAIN = ROOT / "Reference_Data/hg19ToHg38.over.chain.gz"
OUT = ROOT / "Results"
N_PERM = 20_000
SEED = 20261105


def consensus(state: str, positive_fraction: object) -> bool:
    if state not in {"A", "B"} or pd.isna(positive_fraction):
        return False
    value = float(positive_fraction)
    return value >= 0.8 if state == "A" else value <= 0.2


def read_track(path: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(
        path, sep="\t", header=None, names=["chrom", "start", "end", name]
    )


def zscore_by_chrom(frame: pd.DataFrame, column: str) -> pd.Series:
    return frame.groupby("chrom")[column].transform(
        lambda values: (values - values.mean()) / values.std(ddof=1)
    )


def main() -> None:
    primary = read_track(PRIMARY / "GSM4872656_ABcompartment1.bedGraph.gz", "hMSC")
    primary = primary.merge(
        read_track(PRIMARY / "GSM4872658_ABcompartment3.bedGraph.gz", "adipocyte"),
        on=["chrom", "start", "end"], how="inner"
    )
    primary = primary[primary["chrom"].str.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])")].copy()
    primary["hMSC_z"] = zscore_by_chrom(primary, "hMSC")
    primary["adipocyte_z"] = zscore_by_chrom(primary, "adipocyte")
    primary["primary_delta_z"] = primary["adipocyte_z"] - primary["hMSC_z"]
    primary["primary_switch"] = np.sign(primary["hMSC"]) != np.sign(primary["adipocyte"])
    primary["primary_direction"] = np.where(
        primary["primary_switch"],
        np.where(primary["hMSC"] > 0, "A_to_B", "B_to_A"),
        "no_switch",
    )

    lifter = LiftOver(str(CHAIN))
    lifted_rows = []
    for row in primary.itertuples(index=False):
        left = lifter.convert_coordinate(row.chrom, int(row.start))
        right = lifter.convert_coordinate(row.chrom, int(row.end) - 1)
        if len(left) != 1 or len(right) != 1:
            continue
        if left[0][0] != right[0][0] or left[0][2] != "+" or right[0][2] != "+":
            continue
        start38, end38 = int(left[0][1]), int(right[0][1]) + 1
        if end38 <= start38:
            continue
        lifted_rows.append(
            {
                "chrom_hg19": row.chrom,
                "start_hg19": row.start,
                "end_hg19": row.end,
                "chrom": left[0][0],
                "start_grch38": start38,
                "end_grch38": end38,
                "primary_hMSC": row.hMSC,
                "primary_adipocyte": row.adipocyte,
                "primary_delta_z": row.primary_delta_z,
                "primary_switch": bool(row.primary_switch),
                "primary_direction": row.primary_direction,
            }
        )
    lifted = pd.DataFrame(lifted_rows)

    external = pd.read_csv(EXTERNAL)
    external = external[external["chrom"].str.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])")].copy()
    external["tested"] = [
        consensus(state_h, frac_h) and consensus(state_a, frac_a) and pd.notna(delta)
        for state_h, frac_h, state_a, frac_a, delta in zip(
            external["hMSC_state"], external["hMSC_positive_fraction"],
            external["AC_state"], external["AC_positive_fraction"],
            external["AC_vs_hMSC_delta_z"],
        )
    ]
    tested = external[external["tested"]].copy()
    tested["center"] = (tested["start"] + tested["end"]) // 2

    comparison_rows = []
    primary_switch_external_bins: set[tuple[str, int]] = set()
    for row in lifted.itertuples(index=False):
        subset = tested[
            (tested["chrom"] == row.chrom)
            & (tested["center"] >= row.start_grch38)
            & (tested["center"] < row.end_grch38)
        ]
        if subset.empty:
            continue
        external_hmsc = float(subset["hMSC_mean_E1"].mean())
        external_adipocyte = float(subset["AC_mean_E1"].mean())
        external_switch = np.sign(external_hmsc) != np.sign(external_adipocyte)
        external_direction = (
            ("A_to_B" if external_hmsc > 0 else "B_to_A") if external_switch else "no_switch"
        )
        for ext in subset.itertuples(index=False):
            key = (ext.chrom, int(ext.start))
            if row.primary_switch:
                primary_switch_external_bins.add(key)
        comparison_rows.append(
            {
                **row._asdict(),
                "n_external_bins": len(subset),
                "external_hMSC_mean": external_hmsc,
                "external_adipocyte_mean": external_adipocyte,
                "external_delta_z_mean": float(subset["AC_vs_hMSC_delta_z"].mean()),
                "external_switch_aggregated": bool(external_switch),
                "external_direction_aggregated": external_direction,
                "switch_status_agreement": bool(row.primary_switch == external_switch),
                "switch_direction_agreement": bool(
                    row.primary_switch and external_switch and row.primary_direction == external_direction
                ),
            }
        )
    comparison = pd.DataFrame(comparison_rows)
    # Compartment eigenvector signs are arbitrary by chromosome. Orient the
    # primary track against the independent hMSC reference only, then apply the
    # same chromosome sign to both primary lineages. The adipocyte contrast is
    # not used to choose orientation.
    orientation_sign = {}
    for chrom, frame in comparison.groupby("chrom"):
        rho_hmsc = frame["primary_hMSC"].corr(frame["external_hMSC_mean"], method="spearman")
        orientation_sign[chrom] = 1 if pd.isna(rho_hmsc) or rho_hmsc >= 0 else -1
    comparison["primary_orientation_sign"] = comparison["chrom"].map(orientation_sign)
    comparison["primary_hMSC_oriented"] = (
        comparison["primary_hMSC"] * comparison["primary_orientation_sign"]
    )
    comparison["primary_adipocyte_oriented"] = (
        comparison["primary_adipocyte"] * comparison["primary_orientation_sign"]
    )
    comparison["primary_delta_z_oriented"] = (
        comparison["primary_delta_z"] * comparison["primary_orientation_sign"]
    )
    comparison["primary_direction_oriented"] = np.where(
        comparison["primary_switch"],
        np.where(comparison["primary_hMSC_oriented"] > 0, "A_to_B", "B_to_A"),
        "no_switch",
    )
    comparison["switch_direction_agreement"] = (
        comparison["primary_switch"]
        & comparison["external_switch_aggregated"]
        & (comparison["primary_direction_oriented"] == comparison["external_direction_aggregated"])
    )
    comparison.to_csv(OUT / "cross_dataset_1mb_compartment_comparison.csv", index=False)

    rho, rho_p = spearmanr(
        comparison["primary_delta_z_oriented"], comparison["external_delta_z_mean"], nan_policy="omit"
    )
    primary_switches = comparison[comparison["primary_switch"]]
    both_switch = primary_switches[primary_switches["external_switch_aggregated"]]

    rows_by_chrom = {
        chrom: frame.sort_values("start").reset_index(drop=True)
        for chrom, frame in tested.groupby("chrom", sort=False)
    }
    observed = 0
    direction_concordant = 0
    primary_direction_by_external_bin: dict[tuple[str, int], str] = {}
    for row in lifted.itertuples(index=False):
        if not row.primary_switch:
            continue
        sign = orientation_sign.get(row.chrom, 1)
        oriented_hmsc = row.primary_hMSC * sign
        direction = "A_to_B" if oriented_hmsc > 0 else "B_to_A"
        subset = tested[
            (tested["chrom"] == row.chrom)
            & (tested["center"] >= row.start_grch38)
            & (tested["center"] < row.end_grch38)
        ]
        for ext in subset.itertuples(index=False):
            primary_direction_by_external_bin[(ext.chrom, int(ext.start))] = direction
    null = np.zeros(N_PERM, dtype=int)
    rng = np.random.default_rng(SEED)
    for chrom, frame in rows_by_chrom.items():
        selected = frame["AC_vs_hMSC_consensus_effect_difference"].astype(int).to_numpy()
        annotation = np.array(
            [(chrom, int(start)) in primary_switch_external_bins for start in frame["start"]],
            dtype=int,
        )
        observed += int(np.sum(selected * annotation))
        for ext in frame[selected.astype(bool)].itertuples(index=False):
            direction = "A_to_B" if ext.hMSC_state == "A" else "B_to_A"
            if primary_direction_by_external_bin.get((chrom, int(ext.start))) == direction:
                direction_concordant += 1
        counts = np.rint(
            np.fft.ifft(np.fft.fft(selected) * np.conj(np.fft.fft(annotation))).real
        ).astype(int)
        null += counts[rng.integers(0, len(frame), size=N_PERM)]
    empirical_p = (1 + int(np.sum(null >= observed))) / (N_PERM + 1)
    pd.DataFrame({"overlap_count": null}).to_csv(
        OUT / "cross_dataset_overlap_circular_shift_null.csv", index=False
    )
    summary = [
        ("primary_autosomal_1mb_bins", len(primary)),
        ("primary_bins_uniquely_lifted", len(lifted)),
        ("primary_bins_with_external_coverage", len(comparison)),
        ("primary_adipocyte_hMSC_switch_bins", int(primary["primary_switch"].sum())),
        ("mapped_primary_switch_bins", len(primary_switches)),
        ("mapped_primary_switches_also_external_aggregated_switch", len(both_switch)),
        ("mapped_switch_direction_concordant", int(both_switch["switch_direction_agreement"].sum())),
        ("chromosomes_sign_flipped_for_primary_orientation", sum(value == -1 for value in orientation_sign.values())),
        ("delta_z_spearman_rho", float(rho)),
        ("delta_z_spearman_nominal_p", float(rho_p)),
        ("external_selected_bins_overlapping_primary_switch_regions", observed),
        ("external_selected_bins_direction_concordant", direction_concordant),
        ("circular_shift_null_mean_overlap", float(null.mean())),
        ("circular_shift_null_2.5_percentile", float(np.quantile(null, 0.025))),
        ("circular_shift_null_97.5_percentile", float(np.quantile(null, 0.975))),
        ("circular_shift_fold_over_null", float(observed / null.mean())),
        ("circular_shift_empirical_p", empirical_p),
        ("permutations", N_PERM),
        ("seed", SEED),
    ]
    with (OUT / "cross_dataset_reproducibility_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(summary)
    print("\n".join(f"{key}: {value}" for key, value in summary))


if __name__ == "__main__":
    main()
