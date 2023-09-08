import os

# Specify the directory containing your script files
folder_path = "/scratch/devel/pnieto/projects/CSF/analysis/Preprocessing/"

# Loop through each file in the directory
for filename in os.listdir(folder_path):
    if filename.endswith(".Rmd"):  # Adjust the file extension if needed
        file_path = os.path.join(folder_path, filename)
        
        # Read the content of the file
        with open(file_path, "r", encoding="utf-8") as file:
            content = file.read()
        
        # Define the pattern to search for
        old_code = """{r fig.width=9}
Idents(data) <- paste0("RNA_snn_res.", res)
FeaturePlot(
  data,
  pt.size = 1,
  features = "PTPRC",
  cols = c("#ADD8E633", "#E46726"),
  order = TRUE,
  label = TRUE
  )
"""
        
        # Define the replacement code
        new_code = """{r fig.width=9, fig.height = 8}
Idents(data) <- paste0("RNA_snn_res.", res)
FeaturePlot(
  data,
  pt.size = 1,
  features = "PTPRC",
  cols = c("#ADD8E633", "#E46726"),
  order = TRUE,
  label = FALSE
  )
"""
        
        # Replace the old code with the new code
        new_content = content.replace(old_code, new_code)
        
        # Write the modified content back to the file
        with open(file_path, "w", encoding="utf-8") as file:
            file.write(new_content)

print("Code replacement complete.")
