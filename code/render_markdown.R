library(workflowr)
library(rmarkdown)
library(glue)

subprojects <- c("CSF_01", "CSF_02", "CSF_03", "CSF_04", "CSF_05", "CSF_06")
proj_dir <- "/scratch/devel/pnieto/projects/CSF"

# QC
purrr::map(subprojects, function(subproj){
  print(subproj)

  if (subproj == "CSF_06"){
    render(
      "/scratch/devel/pnieto/scripts/01_Overview_and_QC_filtering_template.Rmd",
      output_format = "html_document",
      output_file = glue("{proj_dir}/docs/01_Overview_and_QC_filtering_{subproj}.html"),
      params = list(proj_dir = glue("{proj_dir}/data/{subproj}"), proj_name = subproj,
                    out_dir = glue("{proj_dir}/data/{subproj}/output/01_QC_filtering"),
                    min_lib_size = 500, max_lib_size = 25000, min_n_genes = 200, max_n_genes = 6000, max_pct_mt = 25)
    )
    }
  else {
    render(
      "/scratch/devel/pnieto/scripts/01_Overview_and_QC_filtering_template.Rmd",
      output_format = "html_document",
      output_file = glue("{proj_dir}/docs/01_Overview_and_QC_filtering_{subproj}.html"),
      params = list(proj_dir = glue("{proj_dir}/data/{subproj}"), proj_name = subproj,
                    out_dir = glue("{proj_dir}/data/{subproj}/output/01_QC_filtering"), min_lib_size = 800,
                    max_lib_size = 25000, min_n_genes = 350, max_n_genes = 6000, max_pct_mt = 20)
    )
  }
})

# preliminary clustering
purrr::map(subprojects, function(subproj){
  print(subproj)
  # get patients for each subproject
  patients <- list.dirs(path = glue("{proj_dir}/data/{subproj}/jobs"), full.names = FALSE, recursive = FALSE)
  print(patients)

  purrr::map(patients, function(pat){
    print(pat)
    render(
      glue("{proj_dir}/analysis/Preprocessing/{subproj}_{pat}_Processing_clustering.Rmd"),
      output_format = "html_document",
      output_file = glue("{proj_dir}/docs/{subproj}_{pat}_Processing_clustering.html")
    )
  })
})

# CNV annotation
purrr::map(subprojects, function(subproj){
  print(subproj)
  # get patients for each subproject
  patients <- list.dirs(path = glue("{proj_dir}/data/{subproj}/jobs"), full.names = FALSE, recursive = FALSE)
  print(patients)

  purrr::map(patients, function(pat){
    print(pat)
    render(
      "{proj_dir}/analysis/03_CNV_Analysis_Results.Rmd",
      output_format = "html_document",
      output_file = glue("{proj_dir}/docs/{subproj}_{pat}_CNV_Analysis_Results.html"),
      params = list(subproject = subproj, patient = pat)
    )
  })
})

# high resolution annotation
purrr::map(subprojects, function(subproj){
  print(subproj)
  # get patients for each subproject
  patients <- list.dirs(path = glue("{proj_dir}/data/{subproj}/jobs"), full.names = FALSE, recursive = FALSE)
  print(patients)

  purrr::map(patients, function(pat){
    print(pat)
    render(
      glue("{proj_dir}/analysis/04_Automatic_HR_Annotation.Rmd"),
      output_format = "html_document",
      output_file = glue("{proj_dir}/docs/{subproj}_{pat}_Automatic_HR_Annotation.html"),
      params = list(subproject = subproj, patient = pat)
    )
  })
})

# DLBCL
subprojects <- c("CSF_01", "CSF_01", "CSF_02", "CSF_02", "CSF_03", "CSF_03", "CSF_05", "CSF_05")
patients <- c("4608", "5700", "7921", "7974", "8102", "3054", "8084", "4759")

for (i in 1:length(subprojects)){
  print(i)
  print(patients[i])

  render(
    glue("{proj_dir}/analysis/{subprojects[i]}_{patients[i]}_DLBCL_analysis.Rmd"),
    output_format = "html_document",
    output_file = glue("{proj_dir}/docs/DLBCL_analysis_{subprojects[i]}_{patients[i]}.html")
  )
}

###

# formatting of the webpage
render(glue("{proj_dir}/analysis/00_CSF_QC_cellranger.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/00_CSF_QC_cellranger.html"))

render(glue("{proj_dir}/analysis/index.Rmd"),
    	output_format = "html_document",
    	output_file = glue("{proj_dir}/docs/index.html"))

render(glue("{proj_dir}/analysis/01_filtering.Rmd"),
	output_format = "html_document",
	output_file = glue("{proj_dir}/docs/01_filtering.html"))

render(glue("{proj_dir}/analysis/02_CSF_preliminary_annotation.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/02_CSF_preliminary_annotation.html"))

render(glue("{proj_dir}/analysis/03_CSF_automatic_annotation.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/03_CSF_automatic_annotation.html"))

render(glue("{proj_dir}/analysis/03_CSF_CNV_Results.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/03_CSF_CNV_Results.html"))

render(glue("{proj_dir}/analysis/06_Disease_entity_analysis.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/06_Disease_entity_analysis.html"))

render(glue("{proj_dir}/analysis/Glioblastoma_integration.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/Glioblastoma_integration.html"))

## scVI integration

render(glue("{proj_dir}/analysis/05_CSF_Integration.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/05_CSF_Integration.html"))

render(glue("{proj_dir}/analysis/Bcell_integration.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/Bcell_integration.html"))

render(glue("{proj_dir}/analysis/Tcell_integration.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/Tcell_integration.html"))

render(glue("{proj_dir}/analysis/Myeloid_integration.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/Myeloid_integration.html"))

render(glue("{proj_dir}/analysis/Immunecell_integration.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/Immunecell_integration.html"))


render(glue("{proj_dir}/analysis/B_cell_scVI_integration_analysis.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/B_cell_scVI_integration_analysis.html"))

render(glue("{proj_dir}/analysis/Myeloid_scVI_integration_analysis.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/Myeloid_scVI_integration_analysis.html"))

render(glue("{proj_dir}/analysis/T_cells_scVI_integration_analysis.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/T_cells_scVI_integration_analysis.html"))

render(glue("{proj_dir}/analysis/Immune_cells_scVI_integration_analysis.Rmd"),
       output_format = "html_document",
       output_file = glue("{proj_dir}/docs/Immune_cells_scVI_integration_analysis.html"))


# LISI score
cell_types <- c("B cells", "T cells", "Myeloid cells", "Immune cells")

for (c in cell_types){
  render(
    glue("{proj_dir}/analysis/LISI_score.Rmd"),
    params = list(cell = c),
    output_format = "html_document",
    output_file = glue("{proj_dir}/docs/LISI_score_{c}.html")
  )
}
