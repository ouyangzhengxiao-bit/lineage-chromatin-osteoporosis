# Human lineage-associated chromatin remodeling in osteoporosis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21856116.svg)](https://doi.org/10.5281/zenodo.21856116)

Reproducibility repository for the manuscript **“Human Lineage-Associated Chromatin Compartment Remodeling Prioritizes Osteoporosis Regulators.”**

The study integrates lineage-resolved chromatin compartments with skeletal single-cell transcriptomics, histone-mark profiles, bone-mineral-density genetic annotations, and serum proteomic evidence. The repository contains the analysis and figure-generation code, panel-level source data, supplementary tables, and frozen environment specifications used for the submitted version.

## Main findings represented in this release

- 777 prespecified adipocyte-versus-human-mesenchymal-stromal-cell compartment differences in the independent Hi-C series.
- 12 of 60 curated bone-related loci near these differences, corresponding to a 2.21-fold conditional enrichment over chromosome-wise circular shifts.
- 91 predefined candidate genes, including 76 with at least four measured evidence layers and 11 strict five-layer candidates.
- Coverage-aware prioritization of complementary candidate groups led by TNFRSF11B, COL1A2, SATB2, TWIST1, WLS, SFRP4, SPP1, BGLAP, COL1A1, and BMP4.

## Repository structure

- `Code/` — analysis, audit, source-data export, and figure-generation scripts.
- `Source_Data/` — panel-level source tables used to regenerate the principal displays.
- `Supplementary/` — machine-readable supplementary tables and external Hi-C audit outputs.
- `environment.yml` — frozen Python/Conda environment.
- `R_SESSION_REQUIREMENTS.txt` — R version and package requirements.
- `REPRODUCE.md` — command-level reproduction guide.
- `CITATION.cff` — machine-readable software citation metadata.

Large public input matrices are not duplicated. Retrieve them from the accessions documented in the manuscript: GSE160394, GSE151315, GSE151319/GSE151324, GSE317068, GSE194122, and SRP260624. UK Biobank genetic and proteomic data remain subject to approved-access terms.

## Environment

```bash
conda env create -f environment.yml
conda activate 3d-osteoporosis-reproducibility
```

See [`REPRODUCE.md`](REPRODUCE.md) for the required inputs and commands.

## Citation

Please cite the version-specific DOI for release v1.0.0: [10.5281/zenodo.21856117](https://doi.org/10.5281/zenodo.21856117). The concept DOI [10.5281/zenodo.21856116](https://doi.org/10.5281/zenodo.21856116) always resolves to the latest archived release.

## Authors

- Tingting Tan, MD
- Zhengxiao Ouyang, MD, PhD

## Licensing

The analysis code is distributed under the MIT License. Data files and derived tables are described in [`DATA_LICENSE.md`](DATA_LICENSE.md); underlying third-party datasets remain governed by their original access and reuse terms.
