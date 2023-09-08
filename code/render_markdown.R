library(workflowr)
library(rmarkdown)
library(glue)

subprojects <- c("CSF_01", "CSF_02", "CSF_03", "CSF_04", "CSF_05", "CSF_06")

# QC
purrr::map(subprojects, function(subproj){
  print(subproj)

  if (subproj == "CSF_06"){
    render(
      "/scratch/devel/pnieto/scripts/01_Overview_and_QC_filtering_template.Rmd",
      output_format = "html_document",
      output_file = glue("/scratch/devel/pnieto/projects/CSF/docs/01_Overview_and_QC_filtering_{subproj}.html"),
      params = list(proj_dir = glue("/scratch/devel/pnieto/projects/CSF/data/{subproj}"), proj_name = subproj,
                    out_dir = glue("/scratch/devel/pnieto/projects/CSF/data/{subproj}/output/01_QC_filtering"),
                    min_lib_size = 500, max_lib_size = 25000, min_n_genes = 200, max_n_genes = 6000, max_pct_mt = 25)
    )
    }
  else {
    render(
      "/scratch/devel/pnieto/scripts/01_Overview_and_QC_filtering_template.Rmd",
      output_format = "html_document",
      output_file = glue("/scratch/devel/pnieto/projects/CSF/docs/01_Overview_and_QC_filtering_{subproj}.html"),
      params = list(proj_dir = glue("/scratch/devel/pnieto/projects/CSF/data/{subproj}"), proj_name = subproj,
                    out_dir = glue("/scratch/devel/pnieto/projects/CSF/data/{subproj}/output/01_QC_filtering"), min_lib_size = 800,
                    max_lib_size = 25000, min_n_genes = 350, max_n_genes = 6000, max_pct_mt = 20)
    )
  }
})

# preliminary clustering
purrr::map(subprojects, function(subproj){
  print(subproj)
  # get patients for each subproject
  patients <- list.dirs(path = glue("/scratch/devel/pnieto/projects/CSF/data/{subproj}/jobs"), full.names = FALSE, recursive = FALSE)
  print(patients)

  purrr::map(patients, function(pat){
    print(pat)
    render(
      glue("/scratch/devel/pnieto/projects/CSF/analysis/Preprocessing/{subproj}_{pat}_Processing_clustering.Rmd"),
      output_format = "html_document",
      output_file = glue("/scratch/devel/pnieto/projects/CSF/docs/{subproj}_{pat}_Processing_clustering.html")
    )
  })
})

# CNV annotation
purrr::map(subprojects, function(subproj){
  print(subproj)
  # get patients for each subproject
  patients <- list.dirs(path = glue("/scratch/devel/pnieto/projects/CSF/data/{subproj}/jobs"), full.names = FALSE, recursive = FALSE)
  print(patients)

  purrr::map(patients, function(pat){
    print(pat)
    render(
      "/scratch/devel/pnieto/projects/CSF/analysis/03_CNV_Analysis_Results.Rmd",
      output_format = "html_document",
      output_file = glue("/scratch/devel/pnieto/projects/CSF/docs/{subproj}_{pat}_CNV_Analysis_Results.html"),
      params = list(subproject = subproj, patient = pat)
    )
  })
})

# high resolution annotation
purrr::map(subprojects, function(subproj){
  print(subproj)
  # get patients for each subproject
  patients <- list.dirs(path = glue("/scratch/devel/pnieto/projects/CSF/data/{subproj}/jobs"), full.names = FALSE, recursive = FALSE)
  print(patients)

  purrr::map(patients, function(pat){
    print(pat)
    render(
      glue("/scratch/devel/pnieto/projects/CSF/analysis/04_Automatic_HR_Annotation.Rmd"),
      output_format = "html_document",
      output_file = glue("/scratch/devel/pnieto/projects/CSF/docs/{subproj}_{pat}_Automatic_HR_Annotation.html"),
      params = list(subproject = subproj, patient = pat)
    )
  })
})

###

# formatting of the webpage
render("/scratch/devel/pnieto/projects/CSF/analysis/00_CSF_QC_cellranger.Rmd",
       output_format = "html_document",
       output_file = "/scratch/devel/pnieto/projects/CSF/docs/00_CSF_QC_cellranger.html")

render("/scratch/devel/pnieto/projects/CSF/analysis/index.Rmd",
    	output_format = "html_document",
    	output_file = "/scratch/devel/pnieto/projects/CSF/docs/index.html")

render("/scratch/devel/pnieto/projects/CSF/analysis/01_filtering.Rmd",
	output_format = "html_document",
	output_file = "/scratch/devel/pnieto/projects/CSF/docs/01_filtering.html")

render("/scratch/devel/pnieto/projects/CSF/analysis/02_CSF_preliminary_annotation.Rmd",
       output_format = "html_document",
       output_file = "/scratch/devel/pnieto/projects/CSF/docs/02_CSF_preliminary_annotation.html")

render("/scratch/devel/pnieto/projects/CSF/analysis/03_CSF_automatic_annotation.Rmd",
       output_format = "html_document",
       output_file = "/scratch/devel/pnieto/projects/CSF/docs/03_CSF_automatic_annotation.html")

render("/scratch/devel/pnieto/projects/CSF/analysis/03_CSF_CNV_Results.Rmd",
       output_format = "html_document",
       output_file = "/scratch/devel/pnieto/projects/CSF/docs/03_CSF_CNV_Results.html")

render("/scratch/devel/pnieto/projects/CSF/analysis/05_CSF_Integration.Rmd",
       output_format = "html_document",
       output_file = "/scratch/devel/pnieto/projects/CSF/docs/05_CSF_Integration.html")

dir.create("/scratch/devel/pnieto/projects/CSF/docs/00_CSF_QC_cellranger_files")
dir.create("/scratch/devel/pnieto/projects/CSF/docs/00_CSF_QC_cellranger_files/figure-html")

