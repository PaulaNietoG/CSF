---
title: "B cells scVI integrated"
author:
- name: Paula Nieto
  affiliation: 
  - Centro Nacional de Análisis Genómico (CNAG)
  email: paula.nieto@cnag.crg.eu
date: '`r format(Sys.Date(), "%B %d, %Y")`'
output: 
  html_document:
    toc: true
    toc_float: true
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = FALSE,
  results = "hide",
  warning = FALSE,
  message = FALSE,
  tidy = TRUE
)
```

```{r surce-utils, include=FALSE}
tryCatch({
  source("/scratch/devel/pnieto/scripts/r_utils/utils.R")
}, error = function(e) {
  source("S:/scripts/r_utils/utils.R")
})
library(dittoSeq)
```

```{r set-params}
root <- get_root_dir()
proj_dir <- glue("{root}/CSF")
data_dir <- glue("{proj_dir}/output/integration/B cells")
out_dir <- data_dir
```

```{r create-output-folder}
# create output folder if it doesn't exist yet
if (!file.exists(out_dir)) {
  dir.create(out_dir)
}
```

```{r load and process data}
data <- readRDS(glue("{data_dir}/B_cells_merged.rds"))

data <- data %>%
  NormalizeData() %>% 
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:20, verbose = FALSE, reduction.key = "UMAP_scVI_", reduction.name = "umap_scvi")
```

```{r scVI metadata and UMAP coordinates}
# Read the CSV file with the clusters from scVI
clusters <- read.csv(glue("{data_dir}/B_cells_scVI_clusters.csv"))
# add clusters metadata to seu obj
data <- AddMetaData(data, metadata = clusters)

umap <- read.csv(glue("{data_dir}/B_cells_scVI_UMAP.csv"))
```

```{r}

```

