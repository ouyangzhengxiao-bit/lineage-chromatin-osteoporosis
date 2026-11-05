# Mesenchymal lineage chromatin-compartment reproducibility analysis

Version 2.0.0 is permanently archived at [Zenodo](https://doi.org/10.5281/zenodo.23969522).

This repository release supports the reframed manuscript, “Cross-dataset reproducibility of mesenchymal lineage chromatin compartment remodeling and its relationship to bone mineral density loci.” It replaces the earlier candidate-prioritization interpretation with prespecified, genome-wide reproducibility and overlap tests.

## Main findings reproduced here

- Fifteen retained GSE151319/SRP260624 in situ Hi-C libraries yielded 12,608 analyzable 160-kb bins and 777 hMSC–adipocyte opposite-state, effect-filtered differences.
- Of 361 primary 1-Mb switches with external coverage, 107 were also classified as switches in GSE160394 and 55 had concordant direction.
- Continuous cross-dataset compartment changes had Spearman rho = -0.0701; regional switch overlap was 111 observed versus a circular-shift mean of 96.709 (empirical P = 0.1777).
- Complete Morris et al. eBMD signals were not enriched in the selected compartment-difference bins: 29 direct overlaps versus a circular-shift mean of 44.5701 (empirical P = 0.9941).
- None of 12 prespecified direction-by-process tests survived multiple-testing correction (all FDR at least 0.9965).

These results support a descriptive lineage-associated chromatin-compartment map but do not establish population-level association, causal gene assignment, or clinical prediction.

## Directory structure

- `Code/`: scripts for coordinate mapping, independent Hi-C comparison, circular-shift tests, figures, and manuscript generation.
- `Results/`: complete derived CSV outputs from the reframed analyses.
- `Source_Data/`: the eight supporting-information tables submitted with the manuscript.
- `Figures/`: PDF and PNG versions of Figures 1–4.
- `Reference_Data/`: local-only public inputs downloaded by the user; ignored by Git because the source archives are large or maintained by the original repositories.

## Reproduction order

1. `map_switch_bins_to_genes.py`
2. `compare_independent_hic_datasets.py`
3. `test_complete_morris_loci_overlap.py`
4. `export_prespecified_go_gene_sets.R`
5. `test_prespecified_processes_by_circular_shift.py`
6. `make_reframed_figures.py`

The scripts expect the public source datasets and reference files described in `REPRODUCE.md`. Large public input files are not duplicated in this release.

`Fig3b_replicate_correlations.csv` and `Fig3c_pca_scores.csv` are retained in `Source_Data/` as figure-source tables in addition to the eight numbered supporting-information files.

## Scope change from the earlier release

The 60-gene curated enrichment, 52-gene proximity screen, 91-gene multi-layer score, proteomic and ligand–receptor extensions, hematopoietic multiome comparison, and RSPO3 validation claim are not part of the revised evidence chain. Old outputs should not be cited as results of the reframed manuscript.

## License and citation

Code is released under the MIT License. Public source datasets remain subject to their original terms. Cite the version-specific archive DOI `10.5281/zenodo.23969522` when using this release.
