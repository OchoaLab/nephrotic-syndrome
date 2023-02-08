---
title: |
  | \huge Bristol demographics, replication plan
author: |
  | \large Tiffany Tu, Alejandro Ochoa
#  | \scriptsize ---
#  | \href{http://twitter.com/DrAlexOchoa}{\faTwitter ~DrAlexOchoa}
#  | \href{http://ochoalab.github.io}{\faHome ~ochoalab.github.io}
#  | \href{mailto:alejandro.ochoa@duke.edu}{\faEnvelope ~alejandro.ochoa@duke.edu}
institute: |
  | \normalsize StatGen, Biostatistics & Bioinformatics --- Duke University
date: |
  | \large 2023-02-08 --- NS U01 working group
output: 
  beamer_presentation:
    includes:
      in_header: rmd-beamer-header.tex
    fig_crop: false
    latex_engine: lualatex # for emoji
classoption: "12pt, aspectratio=169"
---

# Overview

- Bristol 
  - All are cases
  - Small sample size (intersect with subtype, ancestry)
    - Only joint (all ancestries) analysis makes sense
  - Have to impute (most candidate loci are missing from raw)
- UK Biobank for controls
  - Overkill considering small number of cases
  - Too big to use it all (n=500,000) with same methods (GMMAT)
  - Can subsample, then do a joint analysis
  - If array genotypes, can impute too
  - WGS is more expensive, perhaps overkill
- GnomAD
  - Retrieve allele counts by ancestry, calculate joint p-value with LRT
  - If it can be automated, could test all suggestive loci this way!

# Bristol demographics

Total $n=$ 590 individuals without filters (age filter further reduces counts).

::: columns
:::: column

| Sex    | Count |    % |
|--------|------:|-----:|
| Male   |   353 | 59.8 |
| Female |   237 | 40.2 |

| Diagnosis       | Count |    % |
|-----------------|------:|-----:|
| SSNS            |   368 | 62.4 |
| SRNS            |   151 | 25.6 |
| NS unclassified |    71 | 12.0 |

::::
:::: column

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   402 | 68.1 |
| Asian          |    84 | 14.2 |
| Unknown        |    73 | 12.4 |
| Black          |    16 |  2.7 |
| Mixed          |    15 |  2.5 |

::::
:::

# Bristol diagnosis subtypes

::: columns
:::: column

SSNS only

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   241 | 65.5 |
| Asian          |    57 | 15.5 |
| Unknown        |    50 | 13.6 |
| Black          |     7 |  1.9 |
| Mixed          |    13 |  3.5 |


::::
:::: column

SRNS only

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   108 | 71.5 |
| Asian          |    17 | 11.3 |
| Unknown        |    17 | 11.3 |
| Black          |     8 |  5.3 |
| Mixed          |     1 |  0.7 |

::::
:::

# Age distribution

\centering
\includegraphics[height=0.75\textheight]{../replication-ukbb/age.pdf}

- 47 individuals missing age
- 265 individuals (44.9%) have age < 18

# Bristol demographics, age < 18 only

Total $n=$ 265 individuals

::: columns
:::: column

| Sex    | Count |    % |
|--------|------:|-----:|
| Male   |   158 | 59.6 |
| Female |   107 | 40.4 |

| Diagnosis       | Count |    % |
|-----------------|------:|-----:|
| SSNS            |   182 | 68.7 |
| SRNS            |    78 | 29.4 |
| NS unclassified |     5 |  1.9 |

::::
:::: column

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   154 | 58.1 |
| Asian          |    49 | 18.5 |
| Unknown        |    42 | 15.8 |
| Black          |    11 |  4.2 |
| Mixed          |     9 |  3.4 |

::::
:::

# Bristol diagnosis subtypes, age < 18 only

::: columns
:::: column

SSNS only

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   103 | 56.6 |
| Asian          |    34 | 18.7 |
| Unknown        |    29 | 15.9 |
| Black          |     7 |  3.8 |
| Mixed          |     9 |  4.9 |

::::
:::: column

SRNS only

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |    48 | 61.5 |
| Asian          |    14 | 17.9 |
| Unknown        |    12 | 15.4 |
| Black          |     4 |  5.1 |
| Mixed          |     0 |  0.0 |

::::
:::

# UKBB costs (3,000 pounds = 3,620.18 USD)

\centering
\includegraphics[height=0.85\textheight]{img/2023-02-08_ukbb-costs.png}

