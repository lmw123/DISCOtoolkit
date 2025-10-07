# DISCOtoolkit 1.1.0 (Jun 11, 2024)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13984774.svg)](https://zenodo.org/badge/DOI/10.5281/zenodo.13984774.svg)

DISCOtoolkit is an R package that provides access to data and tools from the [DISCO database](https://www.immunesinglecell.org/). Its functions include:

- Filter and download DISCO data based on sample metadata and specified cell types
- CELLiD: cell type annotation
- scEnrichment: geneset enrichment using DISCO DEGs

## Requirement

DISCOtoolkit depends on the following packages:

-   [R](https://www.r-project.org/) (\>= 4.0.0)
-   [pbmcapply](https://cran.r-project.org/web/packages/pbmcapply/index.html)
-   [stringr](https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html)
-   [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)
-   [progress](https://cran.r-project.org/web/packages/progress/index.html)


## Installation

``` r
devtools::install_github("git@github.com:JinmiaoChenLab/DISCOtoolkit.git")
```

## Basic Usage

### Filter and download DISCO data

``` r
library(DISCOtoolkit)

# find samples from normal lung tissue and sequenced by 10X Genomics 5' platform
# retain samples containing more than 100 Macrophages(or its children)
metadata = FilterDiscoMetadata(
  sample_id = NULL,
  project_id = NULL,
  tissue = "lung",
  disease = NULL,
  platform = c("10x5'"),
  sample_type = c("control", "adjacent normal"),
  cell_type = "Macrophage", 
  cell_type_confidence = "high", 
  include_cell_type_children = T, 
  min_cell_per_sample = 100
)
### print information ###
# Fetching sample metadata
# Filtering sample
# Fetching cell type information
# Fetching ontology from DISCO database
# 18 samples and 64592 cells were found

# download filtered data into 'disco_data' folder
DownloadDiscoData(metadata, output_dir = "disco_data")
```

### CELLiD

``` r
library(DISCOtoolkit)
library(Seurat)

metadata = FilterDiscoMetadata(
  sample_id = "ERX2757110"
)

DownloadDiscoData(metadata, output_dir = "disco_data")

rna = readRDS("disco_data/ERX2757110.rds")
rna = CreateSeuratObject(rna)
rna = NormalizeData(rna)
rna = FindVariableFeatures(rna)
rna = ScaleData(rna)
rna = RunPCA(rna)
rna = FindNeighbors(rna, dims = 1:10)
rna = FindClusters(rna)

rna_average = AverageExpression(rna)
predict_ct = CELLiDCluster(rna = as.matrix(rna_average$RNA))


# It will download reference data and differential expression gene (DEG) data from DISCO and save them in the 'DISCOtmp' folder by default. You can reuse this data for subsequent CELLiD analyses as follow:

ref_data = readRDS("DISCOtmp/ref_data.rds")
ref_deg = readRDS("DISCOtmp/ref_deg.rds")
predict_ct = CELLiDCluster(rna = as.matrix(rna_average$RNA), ref_data = ref_data, ref_deg = ref_deg)


rna$cell_type = predict_ct$predict_cell_type_1[as.numeric(rna$seurat_clusters)]
rna = RunUMAP(rna, dims = 1:10)
DimPlot(rna, group.by = "cell_type", label = T)
```

### scEnrichment

``` r
markers = FindMarkers(rna, ident.1 = 0, only.pos = T, logfc.threshold = 0.5)
cellid_input = data.frame(gene = rownames(markers), logFC = markers$avg_log2FC)
cellid_res = CELLiDEnrichment(cellid_input)

# also it will download 'ref_geneset.rds' in 'DISCOtmp' folder by default,
# You can reuse this data for subsequent CELLiDEnrichment analyses as follow:
ref = readRDS("DISCOtmp/ref_geneset.rds")
cellid_res = CELLiDEnrichment(cellid_input, reference = ref)
```
