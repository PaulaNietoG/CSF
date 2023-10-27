# CSF project
# read files and select all immune cells for downstream integration with scVI
# create merged object and save counts and metadata

library(workflowr)
library(rmarkdown)
library(glue)
library(magrittr)
library(Seurat)
library(dplyr)

# set path depending if on cluster or local
root_dir <- ifelse(
  grepl(dirname(getwd()), pattern = "/scratch/"),
  "/scratch/devel/pnieto", #true statement
  "S:" # false statement
)
proj_dir <- glue("{root_dir}/projects/CSF")

subprojects <- c("CSF_01", "CSF_02", "CSF_03", "CSF_04", "CSF_05") # , "CSF_06") leaving 6 out because of bad quality/low number of cells

files <- purrr::map(subprojects, function(subproj){
  print(subproj)
  # get patients for each subproject
  patients <- list.dirs(path = glue("{proj_dir}/data/{subproj}/jobs"), full.names = FALSE, recursive = FALSE)
  print(patients)

  dir <- glue("{root_dir}/projects/CSF/data/{subproj}/output/05_CellTypist")

  list.files(dir, full.names = TRUE, recursive = TRUE) %>%
    grep(pattern = "_HR_automatic_annotated_CT.rds", value = TRUE) %>%
    grep(pattern = paste(patients, collapse = "|"), value = TRUE)

}) %>% unlist()

data <- purrr::map(files, function(f){
  tmp <- readRDS(f)
  print(unique(tmp$library))
  print(unique(tmp$annot))
  cell_types <- grep(unique(tmp$annot), pattern = "Myeloid|myeloid|Macro|macro|mono|Mono|DC|T cell|t cell|NK|nk|NKT|B cell|Plasma|plasma", value = TRUE)
  print(cell_types)
  tmp <- subset(tmp, annot %in% cell_types, return.null = TRUE)
  print("***")
  tmp
})
data <- data[!sapply(data,is.null)]
data <- merge(data[[1]], data[2:length(data)])

out_dir <- glue("{proj_dir}/output/integration/Immune cells") %T>%
  dir.create()

# save RDS object
saveRDS(data, glue("{out_dir}/Immune_cells_merged.rds"))
# save counts
write.table(as.matrix(GetAssayData(object = data, slot = "counts")),
            glue("{out_dir}/Immune_cells_counts.csv"),
            sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE)
# save metadata
write.table(
  data@meta.data,
  glue("{out_dir}/Immune_cells_metadata.csv"),
  sep = ',', row.names = TRUE, col.names = TRUE, quote = FALSE
)
