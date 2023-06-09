# This script runs inverCNV (https://github.com/broadinstitute/inferCNV/)

# parameters
subproject <- "CSF_03"
patient <- "3087"
ref_cell_type <- c("T cells 7", "Macrophages 0", "DC 3", "T cells 9", "Macrophages 2")

# Load packages
library(Seurat)
library(tidyverse)
# Sys.setenv(JAGS_HOME="C:/Program Files/JAGS/JAGS-4.3.0") # for windows
library(infercnv)
library(magrittr)
library(glue)
set.seed(1234)

# set path depending if on cluster or local
root_dir <- ifelse(
  grepl(dirname(getwd()), pattern = "/scratch/"),
  "/scratch/devel/pnieto/projects", #true statement
  "S:/projects" # false statement
)
# name of the project directory
proj_dir <- glue::glue("{root_dir}/CSF")

# output path to save data from this step (create dir if it doesn't exist)
out_path <- glue::glue("{proj_dir}/data/{subproject}/output/03_inferCNV/{patient}") %T>%
  dir.create(recursive = TRUE)

# annotated object data path
data_path <- glue("{proj_dir}/data/{subproject}/output/02_processing_clustering/{subproject}_{patient}_annotated.rds")

# directory with inferCNV files
infercnv_data_path <- glue::glue("{proj_dir}/data/inferCNV")

# get date string
date_str <- format(Sys.Date(), "%d_%m_%Y")


path_to_gene_order <- glue::glue("{infercnv_data_path}/gencode_v21_gen_pos.complete.txt")
path_to_save <- glue::glue("{out_path}/{patient}_infercnv_obj.rds")

path_to_cell_annot <- glue::glue("{out_path}/{patient}_cell_annot.txt")
path_to_gene_order_filtered <- glue::glue("{out_path}/{patient}_gene_order_filtered.txt")

cutoff_infercnv <- 0.1 # recommended for 10x data

# Load data
print("Loading data...")
seurat <- readRDS(data_path)
table(seurat$annot_2) %>% as.data.frame() %>% arrange(Freq)
sel <- table(seurat$annot_2) %>% as.data.frame() %>% arrange(Freq) %>% filter(Freq >= 5) %>% pull(Var1) %>% as.character()
# seurat <- subset(seurat, annot_2 %in% sel)

gene_order_file <- read_tsv(
  path_to_gene_order,
  col_names = c("gene", "chromosome", "start", "end")
)

# Write cell annotations file
print("Writing cell annotations file...")
seurat$cell_barcodes <- colnames(seurat)
cell_annotations <- seurat@meta.data[, c("cell_barcodes", "annot_2")]
write_tsv(cell_annotations, path_to_cell_annot, col_names = FALSE)

# Obtain gene ordering file
print("Preparing and saving gene annotation file...")
gene_order_file <- as.data.frame(gene_order_file)
gene_order_file$gene <- str_remove(gene_order_file$gene, "\\|ENSG.*$")
seurat <- subset(
  seurat,
  features = rownames(seurat)[rownames(seurat) %in% gene_order_file$gene]
)
gene_order_file <- gene_order_file[match(rownames(seurat), gene_order_file$gene), ]
write_tsv(gene_order_file, path_to_gene_order_filtered, col_names = FALSE)

# Create inferCNV object
print("Creating inferCNV object...")
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = seurat[["RNA"]]@counts,
  annotations_file = path_to_cell_annot,
  delim = "\t",
  gene_order_file = path_to_gene_order_filtered,
  ref_group_names = ref_cell_type
)

# Run inferCNVs
print("Running inferCNV...")
infercnv_obj <- infercnv::run(
  infercnv_obj,
  # analysis_mode="cells",
  cutoff = cutoff_infercnv,
  out_dir = out_path,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  num_threads = 24,
  HMM = TRUE,
  plot_chr_scale = TRUE,
  BayesMaxPNormal = 0.2
)

# Save
print("Saving...")
saveRDS(infercnv_obj, path_to_save)