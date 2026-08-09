#!/usr/bin/env python3
"""Recompute GC-oriented cis eigenvectors for provenance-complete Hi-C samples."""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import cooler
import cooltools
import numpy as np
import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hic-root", required=True)
    parser.add_argument("--gc-track", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--resolution", type=int, default=160_000)
    parser.add_argument("--min-contacts", type=int, default=10_000_000)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    phasing = pd.read_csv(args.gc_track, sep="\t")
    qc_rows = []

    for mcool in sorted(glob.glob(os.path.join(args.hic_root, "*", "*.mcool"))):
        sample = Path(mcool).stem
        uri = f"{mcool}::resolutions/{args.resolution}"
        clr = cooler.Cooler(uri)
        bins = clr.bins()[:]
        weights = bins["weight"].to_numpy()
        contact_sum = int(clr.info["sum"])
        finite_weight_fraction = float(np.isfinite(weights).mean())
        include = contact_sum >= args.min_contacts and finite_weight_fraction >= 0.75
        group = sample.split("_", 1)[0]
        qc_rows.append(
            {
                "sample": sample,
                "group": group,
                "resolution": args.resolution,
                "contact_sum": contact_sum,
                "nnz": int(clr.info["nnz"]),
                "finite_weight_fraction": finite_weight_fraction,
                "included": include,
                "exclusion_reason": (
                    ""
                    if include
                    else "contact_sum below threshold or insufficient finite balancing weights"
                ),
            }
        )
        if not include:
            print(f"EXCLUDE {sample}: contacts={contact_sum:,}, finite={finite_weight_fraction:.3f}")
            continue

        sample_out = out / f"{sample}.cis.vecs.tsv"
        if sample_out.exists():
            print(f"SKIP {sample}: output exists")
            continue

        # Restrict the phasing track to bins represented by the cooler.
        phase = bins[["chrom", "start", "end"]].merge(
            phasing, on=["chrom", "start", "end"], how="left"
        )
        eigvals, eigvecs = cooltools.eigs_cis(
            clr,
            phasing_track=phase,
            n_eigs=3,
            clr_weight_name="weight",
            sort_metric="spearmanr",
        )
        eigvals.to_csv(out / f"{sample}.cis.lam.tsv", sep="\t", index=False)
        eigvecs.to_csv(sample_out, sep="\t", index=False)
        finite_e1 = int(np.isfinite(eigvecs["E1"]).sum())
        print(f"DONE {sample}: contacts={contact_sum:,}, finite E1={finite_e1:,}")

    pd.DataFrame(qc_rows).to_csv(out / "sample_qc.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
