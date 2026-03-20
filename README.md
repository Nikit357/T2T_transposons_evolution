

# Transposable element–host genome evolutionary arms race revealed by multi-modal epigenomic profiling in a telomere-to-telomere human genome reference

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19052416.svg)](https://zenodo.org/records/19052416)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the computational framework, data processing scripts, and analysis notebooks for the study of the epigenetic dynamics of **3.7 million transposable elements (TEs)** using the complete **T2T-CHM13** human genome assembly.

## Overview

The interaction between TEs and the host genome is a long-standing **evolutionary arms race**. While traditionally difficult to map, the advent of the telomere-to-telomere (T2T) assembly allows for a comprehensive assessment of TE epigenetic states, including previously unresolved centromeric and repetitive regions.

Using the **T2T ENCODE** dataset, this project quantifies the epigenetic impact of TEs across evolutionary time by analyzing **seven epigenomic modalities** across **12 human cell lines**:
* **6 Transposon Classes:** LINE, SINE, LTR, SVA, Helitron, and DNA elements.
* **44 Families and 1,122 Subfamilies**.
* **Epigenomic Modalities:** CTCF binding and histone modifications (H3K4me1, H3K4me3, H3K9ac, H3K27ac, H3K27me3, H3K9me3).

### Key Findings
* **SVA elements** exhibit the strongest signatures of an ongoing arms race, characterized by evasion of H3K9me3-mediated heterochromatinization and increased acquisition of CTCF and enhancer marks.
* Transposon-driven evolution is dominated by **evasion of host heterochromatinization** (H3K9me3 and H3K27me3) and invasion into **CTCF-rich architectural regions**.
* Proliferation history and regulatory potential are linked to TE evolutionary age, demonstrated by **Kimura 2-parameter (K2P)** divergence as a proxy.

## Repository Structure

```text
.
├── scripts/                # Python and bash scripts for data acquisition, mapping and large scale calculations
│
├── plots/                  # svg source plots to build the article figures and the ones not included in the article
│   
├── tables_and_tools/       # csv tables with data and metrics. The main Jupyter notebook pigenetic_modifications_ transposons_before_genes_comapping.ipynb is also here
└── README.md
```

These files can be used to fully reproduce the data analysis done in this paper. Start with data download from UCSC Genome Browser and then run bedtools mapping, and finaly analyze the data using the main notebook.

## Installation & Reproducibility

The analysis was performed in a Jupyter notebook environment using Python 3.

### Dependencies
* **Bioinformatics Tools**: bedtools (v2.29+), SAMtools (v1.10), bigWigToBedGraph.
* **Python Libraries**: pandas, numpy, scipy, seaborn, matplotlib, statsmodels, supervenn.
  
## Data Availability
* **Raw Epigenomic Data**: Enrichment tracks (log ratios of experimental vs. control coverage) were obtained from the T2T ENCODE initiative.
* **TE Annotations**: Coordinates and K2P divergence scores were derived from the T2T RepeatMasker track (CHM13v2.0).

## Methods Summary
* **Evolutionary Age**: K2P CpG-adjusted divergence from the consensus sequence was used as a proxy for evolutionary age.
* **Signal Mapping**: Continuous enrichment signal tracks (BigWig) were intersected with TE genomic intervals using bedtools, and signals were averaged per element.
* **Correlation Analysis**: Spearman’s correlation coefficients were used to detect non-linear dependencies between evolutionary age and epigenetic activity.
* **Significance Testing**: Permutation-based background models (100 iterations) and sigmoid approximations were used to control for stochastic effects in low-copy TE subfamilies.
* **Citation** If you use this code or data in your research, please cite: Nikitin, D. (2026). Transposable element–host genome evolutionary arms race revealed by multi-modal epigenomic profiling in a telomere-to-telomere human genome reference. Manuscript submitted for publication.
## Ethical Statement
This study is purely computational and utilizes publicly available datasets from the ENCODE Consortium. No primary collection of human or animal tissue was involved.

**Author**: Daniil Nikitin

**Affiliation**: Institute of Molecular Biology, National Academy of Science of the Republic of Armenia

**Contact**: danya.nikitin.orel@gmail.com 
