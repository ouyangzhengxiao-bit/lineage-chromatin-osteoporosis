# Reproduction notes

## Public inputs

- GSE151319/SRP260624: primary in situ Hi-C libraries.
- GSE160394: processed 1-Mb A/B compartment bedGraph files (`GSM4872656`, `GSM4872657`, and `GSM4872658`).
- Morris et al. 2019 supplementary workbook: Tables 2 and 6 supply the complete eBMD signals used here.
- GENCODE v49 GRCh38 annotation.
- UCSC `hg19ToHg38.over.chain.gz` for coordinate conversion.
- Gene Ontology gene sets for GO:0001501, GO:0001503, GO:0045444, and GO:0048762.

## Software

Python 3 with pandas, NumPy, SciPy, matplotlib, openpyxl, python-docx, and pyliftover; R with the packages imported by `export_prespecified_go_gene_sets.R` and `run_unbiased_pathway_enrichment.R`.

## Fixed analysis choices

- Primary compartment bins: 160 kb.
- Cross-dataset comparison bins: 1 Mb.
- Circular-shift permutations: 20,000.
- GO tests seed: 20261104.
- eBMD and cross-dataset tests seed: 20261105.
- Empirical upper-tail P value: `(1 + count(null >= observed)) / 20001`.

Circular shifts are restricted to the chromosome-specific analyzable-bin universe. Missing values remain missing and are not treated as zero.

## Release check

Before creating the Zenodo version, run every script against the public inputs, verify that the key values in the repository README match the CSV results, and replace all provisional archive language in the manuscript with the version-specific DOI.
