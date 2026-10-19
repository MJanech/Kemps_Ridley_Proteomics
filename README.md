# Kemp's ridley plasma proteomics analysis

This repository contains the analysis scripts used for the revised manuscript:

**Berardone C, Stacy NI, Owens D, Janech MG. _Proteomic analysis of plasma from rescued immature Kemp's ridley sea turtles (Lepidochelys kempii) exposed to light and heavy degrees of oiling from the Deepwater Horizon oil spill._**  
Journal of Proteome Research, manuscript **pr-2026-00577k**.

The analyses compare lightly oiled turtles (OC2) and heavily oiled turtles (OC4) at admission and after rehabilitation using label-free LC-MS/MS plasma proteomics data.

## Data availability

Raw mass-spectrometry data are available through ProteomeXchange/PRIDE under accession **PXD079061** (DOI: **10.6019/PXD079061**).

The finalized protein-level Supporting Information workbook used by the revised scripts is:

`Supplemental_Tables_1_through_8_Limma_Primary.xlsx`

It contains Tables S1-S8, including the weighted spectral count (wSC) matrix, primary limma-trend results, and nonparametric sensitivity-analysis results. All repository scripts that directly read Tables S1-S8 have been updated to use this revised workbook structure, including the limma, Mann-Whitney, Wilcoxon, PCA/PERMANOVA/PERMDISP, GSEA-ranking, and Figures 8-9 workflows.

The clinical mixed-effects analysis uses a separate finalized workbook:

`Supplemental_Table_S10_Clinical_Data_and_Mixed_Models.xlsx`

because the individual-level clinical data used by that model are provided in Table S10A.

## Protein-level statistical framework

The revised manuscript uses **limma-trend as the primary protein-level differential-abundance analysis**.

Weighted spectral counts are transformed as `log2(wSC + 1)` for limma modeling. The `+1` is a mathematical offset used before logarithmic transformation and is not abundance imputation. Numerical zeros exported from Scaffold are retained as non-detections, and missing protein-abundance values are not imputed.

The Mann-Whitney U analyses for between-oiling-category comparisons and paired Wilcoxon signed-rank analyses for admission-versus-rehabilitation comparisons are retained as **nonparametric sensitivity analyses**. The same comparison-specific protein filters are used for the primary limma-trend analysis and corresponding nonparametric sensitivity analysis.

Benjamini-Hochberg correction is applied separately within each protein-level comparison, with FDR < 0.05 considered significant.

The revised protein-level results are:

| Comparison | Proteins retained | Primary limma-trend FDR < 0.05 | Nonparametric sensitivity FDR < 0.05 |
|---|---:|---:|---:|
| S4: OC2 vs OC4 at admission | 424 | 34 | 17 |
| S5: OC2 vs OC4 after rehabilitation | 328 | 0 | 0 |
| S6: OC2 admission vs rehabilitation | 350 | 145 | 87 |
| S7: OC4 admission vs rehabilitation | 347 | 162 | 153 |

The nonparametric sensitivity findings reproduced by limma-trend were 17/17 for S4, 84/87 for S6, and 141/153 for S7.

The group-by-sampling-occasion limma interaction analysis uses the 306 proteins retained across all four protein-level comparisons. No proteins reached BH FDR < 0.05 for the interaction.

## Analysis environment

The reported limma analysis was generated in **R 4.3.2** using **limma 3.58.1**.

The revised limma workflow was also re-run successfully under **R 4.6.1** with **limma 3.68.5**. The significant-protein counts and conclusions were unchanged. The newer environment produced only very small numerical differences in individual statistics.

The Reactome GSEA reported in the manuscript was generated using:

- R 4.3.2
- Bioconductor 3.18
- `fgsea` 1.28.0
- `msigdbr` 26.1.1
- MSigDB database release `2026.1.Hs`
- `BiocParallel` 1.36.0

The GSEA workflow was re-tested in this environment and exactly reproduced the finalized Supporting Information results. Newer versions of `fgsea` can produce small numerical differences in pathway-level enrichment statistics and adjusted p-values. In validation, `fgsea` 1.38.0 did not reproduce the final pathway calls exactly, whereas `fgsea` 1.28.0 did. The version used for the reported analysis should therefore be used when exact reproduction of Supplemental Table S9 is required.

The clinical mixed-effects models were run in R 4.6.1 using `readxl` 1.5.0, `lme4` 2.0.6, `lmerTest` 3.2.1, and `emmeans` 2.0.4.

## Scripts

- `limma_primary_and_interaction_analysis.R`  
  Primary limma-trend protein-level analyses for Tables S4-S7, standard-limma comparison, nonparametric sensitivity comparison, and group-by-sampling-occasion interaction.

- `MannWhitney_nonparametric_sensitivity_analysis.R`  
  Nonparametric sensitivity analyses for OC2-versus-OC4 comparisons at admission and after rehabilitation (Tables S4-S5).

- `Wilcoxon_paired_nonparametric_sensitivity_analysis.R`  
  Nonparametric paired sensitivity analyses for admission-versus-rehabilitation comparisons within OC2 and OC4 (Tables S6-S7).

- `PCA_PERMANOVA_PERMDISP_Figures4_5.R`  
  PCA, PERMANOVA, PERMDISP, transformation/detection-filter sensitivity analyses, and Figures 4-5.

- `GSEA_nonparametric_ranked_analysis.R`  
  Preranked Reactome GSEA using the comparison-specific nonparametric ranking statistic. This pathway-ranking procedure is separate from the primary limma-trend protein-level differential-abundance analysis.

- `GSEA_collapse_redundant_pathways.R`  
  Reduces significant Reactome pathways to representative main pathways using `fgsea::collapsePathways()`.

- `GSEA_Figures6_7.R`  
  Generates pathway Figures 6 and 7 from the collapsed GSEA results.

- `Figures8_9_violin_plots.R`  
  Generates 600-dpi manuscript Figures 8 and 9 from original, untransformed wSC values and the primary limma-trend BH FDR results in Tables S6-S7.

- `clinical_mixed_models_S2_S10.R`  
  Clinical and morphometric mixed-effects analysis for Tables S2 and S10.

- `modify_fasta_headers.py`  
  Utility used to append gene symbols from an NCBI feature table to matching RefSeq FASTA headers. This utility is included for provenance and is not part of the R statistical workflow.

## PCA, PERMANOVA, and PERMDISP

`PCA_PERMANOVA_PERMDISP_Figures4_5.R` performs the exploratory PCA and Bray-Curtis PERMANOVA/PERMDISP analyses used for Figures 4-5, together with log2(wSC + 1) and detection-filter sensitivity analyses.

Paired admission-versus-rehabilitation permutation tests are restricted within turtle. The four-group analysis is treated as a global/descriptive comparison because each turtle contributes repeated observations.

For the original raw-wSC four-group analysis, the script reproduces:

- PERMANOVA F = 6.861694
- R² = 0.318728
- p = 0.0001
- PERMDISP p = 0.1999

## Reactome GSEA workflow

The pathway workflow is sequential:

```text
GSEA_nonparametric_ranked_analysis.R
    -> GSEA_collapse_redundant_pathways.R
    -> GSEA_Figures6_7.R
```

GSEA uses the complete comparison-specific ranked gene list rather than restricting the analysis to proteins that individually reach FDR < 0.05. The ranking metric is:

```text
rank-biserial correlation × [-log10(original Mann-Whitney/Wilcoxon p-value)]
```

Exact ties in the ranking score are broken deterministically with a Hodges-Lehmann nonparametric effect estimate. No limma statistics are used in the reported GSEA. This does not change the designation of limma-trend as the primary **protein-level differential-abundance** analysis; the nonparametric statistic is used specifically to construct the preranked gene list for pathway analysis.

Using R 4.3.2 and `fgsea` 1.28.0, the workflow reproduces:

| Comparison | Ranked unique genes | Reactome pathways tested | BH FDR < 0.05 | Main pathways after collapse |
|---|---:|---:|---:|---:|
| S4 | 391 | 246 | 7 | 4 |
| S5 | 298 | 165 | 0 | 0 |
| S6 | 319 | 188 | 5 | 3 |
| S7 | 317 | 167 | 21 | 7 |

The pathway collapse uses `fgsea::collapsePathways()` with a conditional enrichment threshold of p = 0.05 and `nperm = 2000`. This is a pathway-dependence criterion and is not a Jaccard-similarity cutoff.

## Running the scripts

Run scripts from the repository working directory. The revised S1-S8 workbook may be supplied explicitly, for example:

```bash
Rscript limma_primary_and_interaction_analysis.R Supplemental_Tables_1_through_8_Limma_Primary.xlsx
```

When no command-line workbook argument is supplied, the revised scripts that directly read Tables S1-S8 default to:

```text
Supplemental_Tables_1_through_8_Limma_Primary.xlsx
```

Example sensitivity-analysis commands are:

```bash
Rscript MannWhitney_nonparametric_sensitivity_analysis.R Supplemental_Tables_1_through_8_Limma_Primary.xlsx
Rscript Wilcoxon_paired_nonparametric_sensitivity_analysis.R Supplemental_Tables_1_through_8_Limma_Primary.xlsx
```

The clinical script instead expects:

```bash
Rscript clinical_mixed_models_S2_S10.R Supplemental_Table_S10_Clinical_Data_and_Mixed_Models.xlsx
```

The GSEA scripts create intermediate result directories that are consumed by the next script in the sequence, so they should be run in the order shown above.

## Repository validation

The revised limma-trend primary analysis, Mann-Whitney sensitivity analysis, paired Wilcoxon sensitivity analysis, PCA/PERMANOVA/PERMDISP, GSEA, pathway-collapse, and Figure 6-9 workflows were re-tested during manuscript revision.

The primary limma-trend workflow reproduced the expected significant-protein counts of **34, 0, 145, and 162** for Tables S4-S7. Figures 8 and 9 use the finalized primary limma-trend BH FDR values and original untransformed wSC plotting values.

The GSEA analysis was specifically re-run under the original R 4.3.2 / Bioconductor 3.18 / `fgsea` 1.28.0 environment to confirm exact agreement with the finalized pathway results.
