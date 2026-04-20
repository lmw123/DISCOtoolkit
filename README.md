# DISCOtoolkit <img src="https://immunesinglecell.com/favicon.ico" align="right" height="40"/>

R client for the [DISCO](https://immunesinglecell.com) single-cell RNA-seq database and analysis platform.

- **22,743 samples · 144,831,560 cells · 40 atlases**
- Covers blood, lung, liver, brain and 30+ other tissues
- Includes healthy controls and 100+ disease conditions

> **Note:** DISCOtoolkit v3 was just released and is under active development — expect frequent updates throughout this month. Some features may be unstable. We apologize for any inconvenience.

## Installation

```r
remotes::install_github("lmw123/DISCOtoolkit")
```

**Dependencies**

| Package | Role |
|---------|------|
| `httr2` | HTTP requests |
| `progress` | Download progress bars |
| `cli` | User-facing messages |

---

## Quick start

```r
library(DISCOtoolkit)
# DISCOtoolkit v3.0.0 — 22,743 samples · 144,831,560 cells · 40 atlases
# last update: 2025-04
# immunesinglecell.com
```

---

## 1. Data access from DISCO

```r
# Filter by tissue and/or disease (returns all matches)
meta <- FilterDiscoMetadata(tissue = "blood", disease = "COVID-19")
meta <- FilterDiscoMetadata(tissue = c("blood", "liver"), limit = 500)
head(meta[, c("sample_id", "tissue", "disease", "cell_number")])

# Filter by sample or project ID
meta <- FilterDiscoMetadata(project_id = "GSE174748")

# What values are available for a filter field?
ListMetadataItem("tissue")
#   value  count
# 1 blood   4231
# 2 lung    1802

# Download H5 files; skips already-downloaded ones
res <- DownloadDiscoData(meta, output_dir = "DISCOtmp")
#   sample_id status                 path
# 1 S01_...      ok  DISCOtmp/S01_....h5
# 2 S02_...  failed                 <NA>

res[res$status == "ok", "path"]   # paths to downloaded files
```

---

## 2. CELLiD — cell type annotation

Annotate clusters using the DISCO CELLiD service (same as the website).

```r
# Server-side scoring (default) — identical to the website
predictions <- CELLiDCluster(seu, atlas = "blood")

# Local scoring — downloads reference once, caches in ~/.disco_cache
predictions <- CELLiDCluster(seu, atlas = "blood", local = TRUE)

#   cluster predicted_cell_type  atlas  score
# 1       0              T cell  blood  0.87
# 2       1              B cell  blood  0.91

cluster_map <- setNames(predictions$predicted_cell_type, predictions$cluster)
seu$cell_type <- setNames(cluster_map[as.character(Seurat::Idents(seu))], colnames(seu))
```

---

## 3. scEnrichment

Tests a gene list against two DISCO references and returns a named list.
Implements the [sc_enrichment](https://www.immunesinglecell.com/tool/sc_enrichment) tool.

```r
# Load built-in example data (same as the website examples)
deg   <- disco_example("deg")    # data.frame: gene + logFC (pancreatic DEGs)
genes <- disco_example("genes")  # character vector (ISG/interferon genes)

res <- CELLiDEnrichment(deg)
head(res)
#        type                                    name overlap n_geneset     pval pval_adj odds_ratio
# 1 Cell Type               Acinar cell                    12        45  1.2e-15  9.4e-13         NA
# 2 Phenotype  type 2 diabetes vs control for ...       8        32  3.1e-06  4.8e-04        11.3

res <- CELLiDEnrichment(genes)
```

---

## Citation

If you use DISCO in your research, please cite:

> **DISCO: a database of Deeply Integrated human Single-Cell Omics data**  
> *Nucleic Acids Research*, 2022.  
> https://doi.org/10.1093/nar/gkab1020

> **Rediscovering publicly available single-cell data with the DISCO platform**  
> *Nucleic Acids Research*, 2025.  
> https://doi.org/10.1093/nar/gkae1108

---

## License

MIT © Mengwei Li
