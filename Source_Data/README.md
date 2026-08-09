# Main-figure source-data map

| Display item | Source-data file |
|---|---|
| Figure 1A–B,D–G | `Fig1_single_cell_embedding_expression.csv.gz` |
| Figure 1C | `Fig1_candidate_gene_dotplot.csv` |
| Figure 2A–D,F–I | `Fig2_replicate_contact_matrices_20kb.csv.gz`; `Fig2_group_median_and_difference_matrices_20kb.csv.gz` |
| Figure 2E,J | `Fig2_compartment_E1_tracks.csv` |
| Figure 2 balance audit | `Fig2_balance_audit.csv` |
| Figure 3A | `Fig3a_HiC_sample_QC.tsv` |
| Figure 3B | `Fig3b_replicate_correlations.csv` |
| Figure 3C | `Fig3c_pca_scores.csv` |
| Figure 3D | `Fig3d_compartment_comparison.csv` |
| Figure 3E | `Fig3e_candidate_locus_comparison.csv` |
| Figure 3F | `Fig3f_gwas_switch_overlap.csv`; `Fig3f_gwas_enrichment.csv` |
| Figure 4A–C | `Fig4_GWAS_3D_proximity.csv` |
| Figure 5A,C | `Fig5_3D_ORS_scores.csv` |
| Figure 5B | `Fig5_measurement_coverage.csv` |
| Figure 5D | `Supplementary_3D_ORS_leave_one_layer_out_gene_ranks.csv` |
| Supplementary Figure S1 | `SuppFig1_proteomic_candidate_intersection.csv` |
| Supplementary Figure S2 | `SuppFig2_ligand_receptor_scores.csv` |
| Supplementary ranking audit | `Supplementary_3D_ORS_leave_one_layer_out_gene_ranks.csv` |

The two `Supplementary_*` count files document donor and anatomical-site composition of the stored GSE317068 object. `Supporting_Audit/` retains cohort-count, cell-level differential-expression, and curated switch-assignment exports used during evidence auditing but not plotted in the revised main Figure 1. Files are UTF-8 text; large per-cell and contact-matrix tables are gzip-compressed CSV.
