# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.0
#   kernelspec:
#     display_name: Python [conda env:.conda-2021-sc-de-benchmark-diffxpy0.6]
#     language: python
#     name: conda-env-.conda-2021-sc-de-benchmark-diffxpy0.6-py
# ---

import os
from nxfvars import nxf
os.environ.setdefault("TF_NUM_THREADS", nxf.task('cpus', "16"))
os.environ.setdefault("TF_LOOP_PARALLEL_ITERATIONS", nxf.task('cpus', "16"))

import scanpy as sc
import diffxpy.api as de
import batchglm
import diffxpy
import pandas as pd

print(batchglm.__version__)
print(diffxpy.__version__)

adata = sc.read_h5ad(nxf.input("adata", "../../data/adata-myeloid.h5ad"))
output_file = nxf.input("output_file", "/dev/null")

# +
# sc.pp.subsample(adata, n_obs=50)
# -

# %%time
# pairwise test does not support covariates, plus it's only a loop anyway.
results = []
groups = adata.obs["leiden"].unique()
for group in groups:
    adata.obs[f"is_leiden_{group}"] = adata.obs["leiden"] == group
    results.append(
        de.test.wald(
            adata,
            formula_loc=f"~ 1 + is_leiden_{group} + dataset + n_genes",
            factor_loc_totest=[f"is_leiden_{group}"],
            as_numeric=["n_genes"],
        )
    )

results

res = pd.concat([
    r.summary().assign(comparison = f"{group} vs. Rest") for group,r in zip(groups, results)
])

res.to_csv(output_file)
