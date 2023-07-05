import scanpy

adata = scanpy.read_h5ad("data/CSF_01/output/02_processing_clustering/CSF_01_4608_annotated.h5ad")
adata

# run ikarus
import urllib.request
import anndata
import pandas as pd
from pathlib import Path
from ikarus import classifier, utils, data

signatures_path = "data/ikarus/signatures.gmt"
model_path = "data/ikarus/core_model.joblib"
out_dir = "data/CSF_01/output/04_ikarus/"

# run model with predefined model and signatures
model = classifier.Ikarus(signatures_gmt=signatures_path, out_dir=out_dir, adapt_signatures=True)
model.load_core_model(model_path)
adata.var['gene_symbol'] = adata.var.index
# adata = data.preprocess_adata(adata)
predictions = model.predict(adata, '4608', save=True)
adata.obs.index.to_series().to_csv('data/CSF_01/output/04_ikarus/4608/index.csv', index=False, header=False)

# make use of cnv information and correct prediction
cnv = pd.read_csv("data/CSF_01/output/03_inferCNV/4608/infercnv.observations.txt", sep=" ")
cnv.columns = cnv.columns.astype("int")
cnv.sort_index(axis=1, inplace=True)
cnv = cnv.T
predictions = model.cnv_correct(adata, '4608', save=True)