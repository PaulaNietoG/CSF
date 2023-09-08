import celltypist
# from celltypist import models
import sys

in_file = sys.argv[1]
out_folder = sys.argv[2]
patient = sys.argv[3]

# #Download all the available models.
# models.download_models()
# #Update all models by re-downloading the latest versions if you think they may be outdated.
# models.download_models(force_update = True)

# Celltyping based on the input of count table
# Predict the identity of each input cell.
predictions = celltypist.annotate(in_file, model = 'Immune_All_Low.pkl', majority_voting = True)

#Export the three results to csv tables.
predictions.to_table(folder = out_folder, prefix = patient)
