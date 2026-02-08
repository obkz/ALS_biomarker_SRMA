# ALS_biomarker_SRMA

This repository contains the R scripts used for the systematic review and meta-analysis (SRMA) of diagnostic and prognostic biomarkers in Amyotrophic Lateral Sclerosis (ALS).

**Publication:**
Obara K, et al. *Eur J Neurol*. 2025 Oct;32(10):e70382. 
"Diagnostic and Prognostic Value of Blood and Cerebrospinal Fluid Biomarkers in Amyotrophic Lateral Sclerosis: A Systematic Review and Meta-Analysis"
**PMID:** [41140053](https://pubmed.ncbi.nlm.nih.gov/41140053/)
**DOI:** [10.1111/ene.70382](https://doi.org/10.1111/ene.70382)

## Structure

- `scripts/`: Analysis scripts and custom functions.
  - `scripts/data_prep/`: Data cleaning and preprocessing pipelines.
- `data/`: Study-level datasets used in the analyses.
  - `diagnosis_prognosis_pub.xlsx`: Master dataset used for all analyses.
  - `SuppleAppendix4.xlsx`: Curated version released as Supplementary Appendix 4 in the published manuscript.
- `output/`: Placeholder for figures and tables.


## Requirements

The analysis was performed using **R (version 4.5.0)**.
The main packages used include: `dplyr`, `ggplot2`, `meta`, etc.
(See each script for specific package usage.)

## Citation

Please cite the original study if you adapt this code or utilize the dataset:
> Obara K, et al. Diagnostic and Prognostic Value of Blood and Cerebrospinal Fluid Biomarkers in Amyotrophic Lateral Sclerosis: A Systematic Review and Meta-Analysis. Eur J Neurol. 2025 Oct 27;32(10):e70382. doi: 10.1111/ene.70382

This is an **Open Access article** distributed under the terms of the [Creative Commons Attribution License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).

## License

- **Code**: The analysis scripts in this repository are licensed under the [MIT License](LICENSE).
- **Data & Text**: The content of the associated publication and the provided datasets are distributed under the [Creative Commons Attribution 4.0 International (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/) license.