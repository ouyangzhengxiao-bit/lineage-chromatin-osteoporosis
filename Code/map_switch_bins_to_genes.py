#!/usr/bin/env python3
"""Map genome-wide adipocyte-versus-hMSC compartment differences to genes.

The analysis is deliberately independent of the historical 52-, 60-, and
91-gene candidate sets. Protein-coding transcription start sites from GENCODE
v49 (GRCh38) define the tested gene universe.
"""

from __future__ import annotations

import csv
import gzip
import re
from bisect import bisect_left
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "Source_Data/S2_Table_primary_160kb_compartment_comparison.csv"
GTF = ROOT / "Reference_Data/gencode.v49.annotation.gtf.gz"
OUT = ROOT / "Results"
OUT.mkdir(parents=True, exist_ok=True)


def attributes(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\S+) "([^"]+)"', text))


def read_genes() -> dict[str, list[dict[str, object]]]:
    genes: dict[str, list[dict[str, object]]] = defaultdict(list)
    with gzip.open(GTF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            info = attributes(fields[8])
            if info.get("gene_type") != "protein_coding":
                continue
            raw_chrom = fields[0]
            chrom = raw_chrom if raw_chrom.startswith("chr") else f"chr{raw_chrom}"
            if chrom not in {f"chr{i}" for i in range(1, 23)}:
                continue
            start, end, strand = int(fields[3]) - 1, int(fields[4]), fields[6]
            tss = start if strand == "+" else end - 1
            genes[chrom].append(
                {
                    "gene_id": info.get("gene_id", "").split(".")[0],
                    "gene_symbol": info.get("gene_name", ""),
                    "tss": tss,
                    "strand": strand,
                }
            )
    for chrom in genes:
        genes[chrom].sort(key=lambda row: int(row["tss"]))
    return genes


def consensus(state: str, positive_fraction: str) -> bool:
    if state not in {"A", "B"} or positive_fraction == "":
        return False
    value = float(positive_fraction)
    return value >= 0.8 if state == "A" else value <= 0.2


def main() -> None:
    genes = read_genes()
    tss = {chrom: [int(row["tss"]) for row in rows] for chrom, rows in genes.items()}
    mapped_rows: list[dict[str, object]] = []
    tested_gene_symbols: set[str] = set()
    selected_gene_symbols: set[str] = set()
    direction_gene_symbols: dict[str, set[str]] = defaultdict(set)
    n_tested = n_selected = 0

    with SOURCE.open(newline="") as handle:
        for row in csv.DictReader(handle):
            chrom = row["chrom"]
            if chrom not in genes:
                continue
            tested = (
                consensus(row["hMSC_state"], row["hMSC_positive_fraction"])
                and consensus(row["AC_state"], row["AC_positive_fraction"])
                and row["AC_vs_hMSC_delta_z"] != ""
            )
            if not tested:
                continue
            n_tested += 1
            selected = row["AC_vs_hMSC_consensus_effect_difference"].lower() == "true"
            if selected:
                n_selected += 1
            start, end = int(row["start"]), int(row["end"])
            center = (start + end) // 2
            positions = tss[chrom]
            left, right = bisect_left(positions, start), bisect_left(positions, end)
            in_bin = genes[chrom][left:right]
            for gene in in_bin:
                tested_gene_symbols.add(str(gene["gene_symbol"]))
                if selected:
                    selected_gene_symbols.add(str(gene["gene_symbol"]))
                    direction = f'{row["hMSC_state"]}_to_{row["AC_state"]}'
                    direction_gene_symbols[direction].add(str(gene["gene_symbol"]))

            idx = bisect_left(positions, center)
            candidates = []
            for j in (idx - 1, idx):
                if 0 <= j < len(positions):
                    candidates.append(genes[chrom][j])
            nearest = min(candidates, key=lambda g: abs(int(g["tss"]) - center))
            direction = (
                f'{row["hMSC_state"]}_to_{row["AC_state"]}' if selected else "not_selected"
            )
            mapped_rows.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "hMSC_state": row["hMSC_state"],
                    "adipocyte_state": row["AC_state"],
                    "delta_z": row["AC_vs_hMSC_delta_z"],
                    "selected_difference": selected,
                    "direction": direction,
                    "genes_with_tss_in_bin": ";".join(
                        sorted({str(g["gene_symbol"]) for g in in_bin})
                    ),
                    "n_genes_with_tss_in_bin": len(in_bin),
                    "nearest_gene": nearest["gene_symbol"],
                    "nearest_gene_id": nearest["gene_id"],
                    "nearest_tss": nearest["tss"],
                    "distance_center_to_nearest_tss": abs(int(nearest["tss"]) - center),
                }
            )

    fields = list(mapped_rows[0])
    with (OUT / "genomewide_compartment_difference_gene_map.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(mapped_rows)

    def write_list(name: str, values: set[str]) -> None:
        with (OUT / name).open("w") as handle:
            for value in sorted(values):
                handle.write(value + "\n")

    write_list("gene_universe_protein_coding_tss_in_tested_bins.txt", tested_gene_symbols)
    write_list("selected_genes_protein_coding_tss_in_difference_bins.txt", selected_gene_symbols)
    for direction, values in direction_gene_symbols.items():
        write_list(f"selected_genes_{direction}.txt", values)

    summary = [
        ("tested_bins_with_consensus_states", n_tested),
        ("selected_opposite_state_effect_bins", n_selected),
        ("protein_coding_genes_in_tested_bins", len(tested_gene_symbols)),
        ("protein_coding_genes_in_selected_bins", len(selected_gene_symbols)),
    ]
    summary.extend(
        (f"protein_coding_genes_{key}", len(value))
        for key, value in sorted(direction_gene_symbols.items())
    )
    with (OUT / "genomewide_mapping_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(summary)
    print("\n".join(f"{key}: {value}" for key, value in summary))


if __name__ == "__main__":
    main()
