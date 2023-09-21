# B_cell_integration_scVI
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi
import torch

import session_info
import warnings

torch.set_float32_matmul_precision('high')

# Setting some parameters
warnings.filterwarnings("ignore")

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

# **Load data**
# We are creating an anndata object from pandas dataframes

# load metadata and counts
metadata = pd.read_csv('/scratch/devel/pnieto/projects/CSF/output/integration/B cells/B_cells_metadata.csv')
counts = pd.read_csv('/scratch/devel/pnieto/projects/CSF/output/integration/B cells/B_cells_counts.csv')
# set the row indices of gene_counts_df to match those of cell_metadata_df using the .reindex method
counts = counts.reindex(columns=metadata.index)
# Create an AnnData object
adata = ad.AnnData(X=np.transpose(counts))
# Assign cell-level metadata
adata.obs = metadata

# convert categorical columns to string (pain in the ass for h5ad files)
# Identify categorical columns
categorical_columns = adata.obs.select_dtypes(include=['category', 'object']).columns
# Convert categorical columns to strings
adata.obs[categorical_columns] = adata.obs[categorical_columns].astype(str)
# let's replace all NaNs with empty spaces because otherwise we get errors *sigh*
adata.obs[categorical_columns] = adata.obs[categorical_columns].fillna("")

# save adata object
adata.write('/scratch/devel/pnieto/projects/CSF/output/integration/B cells/B_cells_merged.h5ad', compression="gzip")

# scVI integration
# Number of latent dimension 
n_latent = 30
scvi.model.SCVI.setup_anndata(adata, batch_key='library')
model = scvi.model.SCVI(adata, n_layers=2, n_latent=n_latent, gene_likelihood="nb")
model.train(early_stopping = True)
adata.obsm['X_scVI'] = model.get_latent_representation()

# Save the results
adata.write('/scratch/devel/pnieto/projects/CSF/output/integration/B cells/B_cells_scVI_integrated.h5ad', compression="gzip")
model.save('/scratch/devel/pnieto/projects/CSF/output/integration/B cells/model')