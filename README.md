# Kemp's ridley plasma proteomics analysis

This repository contains the analysis scripts used for the manuscript:

**Berardone C, Stacy NI, Owens D, Janech MG. _Proteomic analysis of plasma from rescued immature Kemp's ridley sea turtles (Lepidochelys kempii) exposed to light and heavy degrees of oiling from the Deepwater Horizon oil spill._**  
Journal of Proteome Research, manuscript **pr-2026-00577k**.

The analyses compare lightly oiled turtles (OC2) and heavily oiled turtles (OC4) at admission and after rehabilitation using label-free LC-MS/MS plasma proteomics data.

## Data availability

Raw mass-spectrometry data are available through ProteomeXchange/PRIDE under accession **PXD079061** (DOI: **10.6019/PXD079061**).

The finalized protein-level Supporting Information workbook used by the scripts in this repository is:

`Supplemental_Tables_1_through_8.xlsx`

It contains Tables S1-S8, including the weighted spectral count (wSC) matrix and the primary protein-level statistical results.

The clinical mixed-effects analysis uses a separate finalized workbook:

`Supplemental_Table_S10_Clinical_Data_and_Mixed_Models.xlsx`

because the individual-level clinical data used by that model are provided in Table S10A.

## Analysis environment

All statistical and proteomic analyses were performed in **R 4.3.2**, except for the clinical mixed-effects models, which were performed in **R 4.6.1**.

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

- `MannWhitney_protein_level_analysis.R`  
  Primary OC2 vs OC4 protein-level comparisons at admission and after rehabilitation (Tables S4-S5).

- `Wilcoxon_paired_protein_level_analysis.R`  
  Primary paired admission vs rehabilitation protein-level comparisons within OC2 and OC4 (Tables S6-S7).

- `limma_sensitivity_and_interaction_analysis.R`  
  limma sensitivity analyses and group-by-sampling-occasion interaction analysis. The manuscript sensitivity analysis used `limma` 3.58.1.

- `PCA_PERMANOVA_PERMDISP_Figures4_5.R`  
  PCA, PERMANOVA, PERMDISP, transformation/detection-filter sensitivity analyses, and Figures 4-5.

- `GSEA_nonparametric_ranked_analysis.R`  
  Nonparametric preranked Reactome GSEA based on the Mann-Whitney U and paired Wilcoxon results.

- `GSEA_collapse_redundant_pathways.R`  
  Reduces significant Reactome pathways to representative main pathways using `fgsea::collapsePathways()`.

- `GSEA_Figures6_7.R`  
  Generates revised pathway Figures 6 and 7 from the collapsed GSEA results.

- `Figures8_9_violin_plots.R`  
  Generates Figures 8 and 9 using original, untransformed wSC values and the BH-adjusted paired Wilcoxon results from Tables S6-S7.

- `clinical_mixed_models_S2_S10.R`  
  Clinical and morphometric mixed-effects analysis for Tables S2 and S10.

- `modify_fasta_headers.py`  
  Utility used to append gene symbols from an NCBI feature table to matching RefSeq FASTA headers. This utility is included for provenance and is not part of the R statistical workflow. A version is also available in the `addGeneSymbol2RefSeq` directory of the `proteomicScripts` repository: https://github.com/neely/proteomicScripts/tree/main/addGeneSymbol2RefSeq

## Primary statistical workflow

Protein wSC values were analyzed using the original nonparametric framework reported in the manuscript. Between-oiling-category comparisons used two-sided Mann-Whitney U tests, and paired admission-versus-rehabilitation comparisons used paired Wilcoxon signed-rank tests. Benjamini-Hochberg correction was applied separately within each comparison.

Numerical zeros in the wSC matrix were retained as Scaffold-exported non-detections. Missing protein-abundance values were not imputed. Comparison-specific sparsity filters are defined in the corresponding analysis scripts.

The primary analysis reproduces the following checkpoints:

| Comparison | Proteins retained | BH FDR < 0.05 |
|---|---:|---:|
| S4: OC2 vs OC4 at admission | 424 | 17 |
| S5: OC2 vs OC4 after rehabilitation | 328 | 0 |
| S6: OC2 admission vs rehabilitation | 350 | 87 |
| S7: OC4 admission vs rehabilitation | 347 | 153 |

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

GSEA uses the complete comparison-specific ranked gene list rather than restricting the analysis to proteins that individually reach FDR < 0.05. The primary ranking metric is:

```text
rank-biserial correlation × [-log10(original Mann-Whitney/Wilcoxon p-value)]
```

Exact ties in the primary score are broken deterministically with a Hodges-Lehmann nonparametric effect estimate. No limma statistics are used in the reported GSEA.

Using R 4.3.2 and `fgsea` 1.28.0, the workflow reproduces:

| Comparison | Ranked unique genes | Reactome pathways tested | BH FDR < 0.05 | Main pathways after collapse |
|---|---:|---:|---:|---:|
| S4 | 391 | 246 | 7 | 4 |
| S5 | 298 | 165 | 0 | 0 |
| S6 | 319 | 188 | 5 | 3 |
| S7 | 317 | 167 | 21 | 7 |

The pathway collapse uses `fgsea::collapsePathways()` with a conditional enrichment threshold of p = 0.05 and `nperm = 2000`. This is a pathway-dependence criterion and is not a Jaccard-similarity cutoff.

## Running the scripts

Run scripts from the repository working directory. For scripts that accept a workbook argument, the finalized S1-S8 workbook may be supplied explicitly, for example:

```bash
Rscript MannWhitney_protein_level_analysis.R Supplemental_Tables_1_through_8.xlsx
```

When no command-line workbook argument is supplied, the S1-S8 scripts default to:

```text
Supplemental_Tables_1_through_8.xlsx
```

The clinical script instead expects:

```bash
Rscript clinical_mixed_models_S2_S10.R Supplemental_Table_S10_Clinical_Data_and_Mixed_Models.xlsx
```

The GSEA scripts create intermediate result directories that are consumed by the next script in the sequence, so they should be run in the order shown above.

## Repository validation

The R scripts using Tables S1-S8 were re-tested during preparation of the manuscript. The Mann-Whitney, paired Wilcoxon, limma, PCA/PERMANOVA/PERMDISP, GSEA, pathway-collapse, and Figure 6-9 workflows completed successfully and reproduced their expected headline results. The GSEA analysis was specifically re-run under the original R 4.3.2 / Bioconductor 3.18 / `fgsea` 1.28.0 environment to confirm exact agreement with the finalized pathway results.
