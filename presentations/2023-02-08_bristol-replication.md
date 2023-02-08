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
| SSNS            |   350 | 59.3 |
| SRNS            |   172 | 29.2 |
| NS unclassified |    68 | 11.5 |

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
| White          |   229 | 65.4 |
| Asian          |    57 | 16.3 |
| Unknown        |    44 | 12.6 |
| Black          |     7 |  2.0 |
| Mixed          |    13 |  3.7 |


::::
:::: column

SRNS only

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   122 | 70.9 |
| Asian          |    23 | 13.4 |
| Unknown        |    18 | 10.5 |
| Black          |     8 |  4.7 |
| Mixed          |     1 |  0.6 |

::::
:::

# Age distribution

\centering
\includegraphics[height=0.75\textheight]{../replication-bristol/age.pdf}

- 47 individuals missing age
- 294 individuals (49.8%) have age < 22

# Bristol demographics, age < 22 only

Total $n=$ 294 individuals

::: columns
:::: column

| Sex    | Count |    % |
|--------|------:|-----:|
| Male   |   178 | 60.5 |
| Female |   116 | 39.5 |

| Diagnosis       | Count |    % |
|-----------------|------:|-----:|
| SSNS            |   191 | 65.0 |
| SRNS            |    96 | 32.7 |
| NS unclassified |     7 |  2.4 |

::::
:::: column

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   174 | 58.8 |
| Asian          |    56 | 19.0 |
| Unknown        |    44 | 15.0 |
| Black          |    11 |  3.7 |
| Mixed          |    10 |  3.4 |

::::
:::

# Bristol diagnosis subtypes, age < 22 only

::: columns
:::: column

SSNS only

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |   110 | 57.6 |
| Asian          |    40 | 20.9 |
| Unknown        |    25 | 13.1 |
| Black          |     7 |  3.7 |
| Mixed          |     9 |  4.7 |

::::
:::: column

SRNS only

| Race/Ethnicity | Count |    % |
|----------------|------:|-----:|
| White          |    59 | 61.5 |
| Asian          |    18 | 18.8 |
| Unknown        |    14 | 14.6 |
| Black          |     4 |  4.2 |
| Mixed          |     1 |  1.0 |

::::
:::

# UKBB costs (3,000 pounds = 3,620.18 USD)

\centering
\includegraphics[height=0.85\textheight]{img/2023-02-08_ukbb-costs.png}

