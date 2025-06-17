# 🔬 CRISPR-Driven Discovery of TMZ Sensitivity Mechanisms in Glioblastoma

This project analyzes genome-wide CRISPR knockout screens across glioblastoma cell lines to identify genes involved in sensitivity or resistance to **Temozolomide (TMZ)**. Using functional enrichment, pathway profiling, and co-dependency analysis, we uncover key biological processes and gene targets relevant to drug response.

---

## Contents

- `scripts/`: R scripts for loading, preprocessing, and analysis
- `data/`: Example CRISPR screen data (e.g., `mmc4.xlsx`)
- `figures/`: Plots (volcano, heatmap, GSEA, etc.)
- `results/`: CSV files for hits, enrichment results, etc.

---

## Dataset

- CRISPR knockout Z-scores for 3 GBM cell lines: **G361, G510, G523**
- MGMT methylation scores
- Gene annotation via Ensembl
- Pathway data: GO, KEGG, Reactome, MSigDB

---

## Project Objectives

1. Identify **sensitizer** and **resistant** genes based on ΔZ scores.
2. Perform **GO, KEGG, Reactome**, and **GSEA** enrichment.
3. Integrate **MGMT methylation** with CRISPR results.
4. Highlight known **DNA repair genes** and **druggable targets**.
5. Visualize key results with **volcano plots**, **heatmaps**, and **ridgeplots**.

---

## 🚀 Getting Started

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/crispr-tmz-gbm.git
cd crispr-tmz-gbm
