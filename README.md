## Citation

> Tu, T., Ochoa, A., Sood, A., Drabik, A., Chryst-Stangl, M., Lane, B., Wu, G., Donovan, F., Harper, U., Chandrasekharappa, S., Esezobor, C., Solarin, A., Hooper, D., Sethna, C., Amaral, S., Kallash, M., Rheault, M., Verghese, P., Dharnidharka, V., Salmon, E., Weng, P., Srivastava, T., Seifert, M.E., Pruette, C., Selewski, D., Gibson, K., Hunley, T., Abeyagunawardena, A., Thalgahagoda, S., Bagga, A., Sinha, A., Webb, N., Greenbaum, L., Gharavi, A., Kiryluk, K., Kretzler, M., Guay-Woodford, L., Sanna-Cherchi, S., Bierzynska, A., Koziell, A., Welsh, G., Saleem, M., Rotimi, C., Chambers, E., Chan, C., CureGN Consortium, PNRC Glomerular disease group, CIBMTR/NMDP Consortium, Jackson, A., Adeyemo, A., Gbadegesin, R., 2025. [Polygenic Risk Scores and HLA Class II Variants are Biomarkers of Corticosteroid Response in Childhood Nephrotic Syndrome](https://www.kidney-international.org/article/S0085-2538(26)00134-1/fulltext). Kidney International. https://doi.org/10.1016/j.kint.2026.01.026

## Overview
1. Genotype quality control and ancestry inference on the discovery cohort.
2. Genome-wide and gene-based association testing for NS, SSNS, and SRNS.
3. SNP annotation and GWAS Catalog submission preparation.
4. Replication in two independent cohorts (Bristol/UKBB and CureGN), including
   ancestry-stratified analyses (AFR, EUR, SAS).
5. Polygenic risk score (PRS) 

Trait and ancestry labels used throughout the code follow the conventions in
`gwas-catalog/README`:

| Label | Meaning |
| --- | --- |
| `NS` | nephrotic syndrome |
| `SSNS` | steroid-sensitive NS |
| `SRNS` | steroid-resistant NS |
| suffix `_AFR` / `_EUR` / `_SAS` | ancestry-stratified subanalysis |
| no suffix | combined multi-ancestry analysis |

## Repository layout

```
quality_control/       Genotype QC, kinship, missingness, admixture, 1KGP merge
association_analyses/  GWAS (SAIGE), SKAT, conditional and meta-analyses,
                       top-SNP annotation
snp_annotation/        Variant annotation, rsID assignment, Excel export
gwas-catalog/          Files and scripts used to prepare the GWAS Catalog
                       submission
replication-bristol/   Replication in the Bristol cohort using UKBB controls
replication-curegn/    Replication in the CureGN cohort using gnomAD controls
prs/                   Polygenic risk score pipeline: LDpred2 (inf, grid, auto),
                       lassosum, SAIGE sumstats prep, PRS evaluation and plots
scripts/               Shared utilities: phenotype/demographic prep, power
                       calculations, BIM intersection, single-cell onset-vs-
                       remission analysis
```
