# READ ME

<a href="https://creativecommons.org">Untitled</a> © 1999 by <a href="https://creativecommons.org">Jane Doe</a> is licensed under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY 4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">

# Drug–Target Interaction Landscape by Gene Essentiality and Disease Ontology

This repository presents a comprehensive analysis of drug–gene interactions by integrating clinical drug data, gene essentiality classifications (FUSIL), and disease ontologies.

---

# Please RUN the ontology mapping first given that some scripts depends on the output.

## Data Sources

| Source            | Description                                     |
|------------------|-------------------------------------------------|
| **Open Targets** | Drug–target–indication associations             |
| **DrugBank**      | Curated drug–gene interaction dataset           |
| **FUSIL**         | Functional gene classification by essentiality |
| **EFO Ontology**  | Disease classes and hierarchical mappings      |

---

## Key Features

- **Comparison of Approved vs Withdrawn Drugs**  
  Evaluate drug distributions by status (e.g. withdrawn, completed)

- **Network-Level Analysis**  
  Degree distributions: drugs per target and targets per drug

- **FUSIL Stratification**  
  Analyze how essentiality bins (CL, DL, SV, VP, VnP) relate to drug targeting

- **Indication-Weighted Targeting Metrics**  
  Calculate indication × drug impact and drug–indication ratios per gene

- **Disease Class Enrichment**  
  Map drug targets to EFO disease classes and visualize FUSIL-based distributions

---

## Visualizations

- Bar plots of unique drugs per dataset  
- Log-log degree distributions of bipartite drug–gene networks  
- Density plots of drugs per gene  
- Indication-weighted drug impact by FUSIL bin  
- Faceted bar charts of disease class enrichment per FUSIL category  

> All plots are generated using `ggplot2`.

---

## 📦 Dependencies

Make sure you have the following R packages installed:

```r
install.packages(c("tidyverse", "readr", "ggplot2", "dplyr", "forcats", "tidyr"))
