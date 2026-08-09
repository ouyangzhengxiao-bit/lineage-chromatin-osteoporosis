# Reproducibility instructions

The submission package contains panel-level source data but does not duplicate
the large public h5ad and multiresolution Hi-C matrices. Download or generate
those inputs from the accessions listed in the manuscript.

## Environment

Create the Python environment from `environment.yml`. The R figure script was
audited with the versions in `R_SESSION_REQUIREMENTS.txt`.

## Required input locations

Set paths explicitly; no script depends on the original author's home or
external-drive layout.

```bash
export EXTERNAL_HIC_WORK=/path/to/work
export GSE317068_H5AD=/path/to/gse317068_processed.h5ad
export HIC_MCOOL_ROOT=/path/to/hic_output
export HIC_E1_GROUP_TABLE=/path/to/group_compartment_comparison.csv
export HIC_GC_TRACK=/path/to/hg38_160kb_gc_track.tsv
```

`HIC_MCOOL_ROOT` must contain `SAMPLE/SAMPLE.mcool` for hMSC_rep1–7,
AC_rep1–6, OB_rep1, and OB_rep6.

## Figure and audit commands

```bash
python Code/run_external_hic_eigs.py \
  --hic-root "$HIC_MCOOL_ROOT" \
  --gc-track "$HIC_GC_TRACK" \
  --out "$EXTERNAL_HIC_WORK/external_validation_160kb/eigs"
python Code/analyze_external_hic.py
python Code/make_singlecell_figure.py
python Code/make_hic_contact_figure.py
python Code/make_external_validation_figure.py
python Code/rebuild_3d_ors_coverage_aware.py \
  --scores Supplementary/Supplementary_Table_2_3D_ORS_scores.csv \
  --h5ad "$GSE317068_H5AD"
python Code/audit_3d_ors_robustness.py \
  --scores Supplementary/Supplementary_Table_2_3D_ORS_scores.csv
Rscript Code/make_core_figures.R
python Code/make_graphical_abstract.py
python Code/freeze_curated_bone_locus_set.py \
  --loci "$EXTERNAL_HIC_WORK/external_validation_160kb/analysis/curated_gwas_loci_external_switch_overlap.csv" \
  --comparison "$EXTERNAL_HIC_WORK/external_validation_160kb/analysis/group_compartment_comparison.csv"
```

The packaged Figure 3 and core-figure scripts read only `Source_Data`. Raw Hi-C
reprocessing commands above are optional audit steps for regenerating those
source tables. The Hi-C locus script requires balanced 20-kb matrices and terminates if
balanced retrieval fails. It does not silently substitute unbalanced contacts.
