
# script to annotate automatically hight res clusters

library(tidyverse)
library(glue)
library(magrittr)
library(Seurat)

subproject = "CSF_01" # args[1]
patient = "4608" # args[2]

# set path depending if on cluster or local
root_dir <- ifelse(
  grepl(dirname(getwd()), pattern = "/scratch/"),
  "/scratch/devel/pnieto", #true statement
  "S:" # false statement
)

# name of the project directory
proj_dir <- glue::glue("{root_dir}/projects/CSF")
# name of output folder
out_dir <- glue("{proj_dir}/data/{subproject}/output/04_high_res_annotation") %T>%
  dir.create()

# load cell type/list
source(glue("{root_dir}/scripts/r_utils/marker_genes.R"))

data <- readRDS(glue("{proj_dir}/data/{subproject}/output/02_processing_clustering/{subproject}_{patient}_annotated.rds"))

# compute signature scores con UCell
data <- UCell::AddModuleScore_UCell(data, features = marker_genes, ncores = 4, name = "")

# get average by "annot_2" cluster
avg <- data@meta.data[, c("annot_2", make.names(names(marker_genes)))]
rownames(avg) <- NULL
averages <- aggregate(. ~ annot_2, data = avg, FUN = mean)
df_avg <- pivot_longer(averages, cols = colnames(averages)[colnames(averages) != "annot_2"])

# heatmap
hm <- ggplot(data = df_avg, aes(x = annot_2, y = name, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +  # Set the color scale from white to red
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for better readability
ggsave(plot = hm, file = glue("{out_dir}/{subproject}_{patient}_heatmap.png"), width = 12)


## assign to each cluster, the cell type with the highest jaccard index
# Initialize an empty list to store the results
highest_column_list <- list()
# Loop through the rows (clusters) of the 'jaccard_indices' matrix
for (row_name in rownames(jaccard_indices)) {
  # Get the row values for the current cluster
  row_values <- jaccard_indices[row_name, ]

  # Check if all values in the row are 0
  if (all(row_values == 0)) {
    highest_column <- "unknown"
  } else {
    # Find the column names with the highest value(s) in the row
    highest_columns <- names(row_values)[row_values == max(row_values)]

    # Check if there is a tie (more than one column with the same highest value)
    if (length(highest_columns) > 1) {
      highest_column <- "unknown"
    } else {
      # Assign the highest_column to the current cluster in the list
      highest_column <- highest_columns
    }
  }

  # Assign the highest_column to the current cluster in the list
  highest_column_list[[row_name]] <- highest_column
}

# assign annotation to annot 2
Idents(data) <- "annot_2"
data <- RenameIdents(data, highest_column_list)
data$auto_annot <- Idents(data)
data$auto_annot <- as.character(data$auto_annot)
# fix tumor clusters if any
data$auto_annot[data$annot %in% grep(unique(data$annot), patern = "Tumor", value = TRUE)] <- "Tumor cells"

# save annotated object
saveRDS(data, glue("{out_dir}/{subproject}_{patient}_HR_automatic_annotated.rds"))

p1 <- DimPlot(data, group.by = c("annot_2", "auto_annot"), pt.size = 1.5, cols = as.vector(pals::polychrome())) +
  labs(
    title = "High Resolution Annotation",
    subtitle = "automatic cell type assignment",
    caption = glue("{subproject} {patient}")
  )
ggsave(plot = p1, file = glue("{out_dir}/{subproject}_{patient}_umaps.png"), width = 12)
