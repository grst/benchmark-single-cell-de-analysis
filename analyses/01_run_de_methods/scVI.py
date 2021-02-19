# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.0
#   kernelspec:
#     display_name: SSH apollo-15 scVI GPU (ssh)
#     language: ''
#     name: rik_ssh_apollo_15_scvigpussh
# ---

import scanpy as sc
import scvi
from nxfvars import nxf

adata = sc.read_h5ad(nxf.input("adata", "../../data/adata-myeloid.h5ad"))
output_file = nxf.input("output_file", "/dev/null")

scvi.data.setup_anndata(adata, batch_key="dataset")

model = scvi.model.SCVI(adata,
                        n_hidden=128,
                        n_layers=2,
                        gene_likelihood='nb',
                        dispersion='gene-batch'
                        )

model.train()

res = model.differential_expression(adata, groupby="leiden")

res

res.sort_values("proba_de", ascending=False).query("comparison == '0 vs Rest'").head(20)

res.to_csv(output_file)
