# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.0
#   kernelspec:
#     display_name: Python [conda env:.conda-2021-sc-de-benchmark-edger]
#     language: python
#     name: conda-env-.conda-2021-sc-de-benchmark-edger-py
# ---

from nxfvars import nxf
cpus = nxf.task("cpus", "16")

# %env MKL_NUM_THREADS={cpus}
# %env NUMEXPR_NUM_THREADS={cpus}
# %env OMP_NUM_THREADS={cpus}
# %env OMP_NUM_cpus={cpus}
# %env MKL_NUM_cpus={cpus}
# %env OPENBLAS_NUM_THREADS={cpus}
# %env OPENBLAS_NUM_cpus={cpus}

# +
import scanpy as sc
import rpy2.rinterface_lib.callbacks
import logging

# from rpy2.robjects import pandas2ri
import anndata2ri
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# +
# Ignore R warning messages
# Note: this can be commented out to get more verbose R output
# rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Automatically convert rpy2 outputs to pandas dataframes
anndata2ri.activate()
# %load_ext rpy2.ipython

plt.rcParams["figure.figsize"] = (8, 8)  # rescale figures
sc.settings.verbosity = 3
# sc.set_figure_params(dpi=200, dpi_save=300)b
sc.set_figure_params(figsize=(5, 5))

# + magic_args="-i cpus" language="R"
# library(dplyr)
# library(edgeR)
# library(BiocParallel)
#
# options(mc.cores=cpus)
# register(MulticoreParam(cpus))
# print(registered())
# -

adata = sc.read_h5ad(nxf.input("adata", "../../data/adata-myeloid.h5ad"))
output_file = nxf.input("output_file", "/dev/null")

# +
# sc.pp.subsample(adata, n_obs=50)
# -

obs = adata.obs
gene_symbols = adata.var_names
counts = adata.X.T.toarray()

# + language="R"
# #' make contrasts: one against all others.
# #' 
# #' @param design design matrix
# #' @param col_data colData or the SingleCellExperiment object. 
# #' @param column column name that is used for the contrasts. Also needs to be 
# #'    specified as first variable in the model.matrix.
# make_contrasts = function(design, col_data, column) {
#     n_clus = length(unique(col_data[[column]]))
#     upper_block = matrix(rep(-1/(n_clus-1), n_clus^2), ncol=n_clus)
#     diag(upper_block) = rep(1, n_clus)
#     lower_block = matrix(rep(0, (ncol(design)-n_clus) * n_clus), ncol=n_clus)
#     contrasts = rbind(upper_block, lower_block)
#     rownames(contrasts) = colnames(design)
#     colnames(contrasts) = colnames(design)[1:n_clus]
#     contrasts
# }

# + magic_args="-i obs -i counts -i gene_symbols" language="R"
# var_data = data.frame(gene_symbol=gene_symbols)
# dge = DGEList(counts=counts, samples=obs, genes=var_data)
# design = model.matrix(~ 0 + leiden + dataset + n_genes, data=obs)
# contrasts = make_contrasts(design, obs, "leiden")

# + language="R"
# head(design)

# + language="R"
# contrasts

# + language="R"
# message("Calculating NormFactors...")
# dge <- calcNormFactors(dge)
#
# message("Estimating Dispersion...")
# dge <- estimateDisp(dge, design = design)
#
# message("Fitting linear model...")
# fit <- glmQLFit(dge, design = design)
#
# message("Testing contrasts...")
# qlfs = sapply(colnames(contrasts), function(cluster) {
#     message(paste0("working on ", cluster))
#     glmQLFTest(fit, contrast=contrasts[,cluster])
# }, USE.NAMES=TRUE, simplify=FALSE)
#
#
# message("Preparing results...")
# tts = sapply(qlfs, function(qlf) {
#     topTags(qlf, n=Inf, adjust.method="BH")
# }, USE.NAMES=TRUE, simplify=FALSE)
#
# all_results = bind_rows(lapply(names(tts), function(cluster) {
#     tts[[cluster]]$table %>% mutate(cluster=cluster)
# })) %>% as_tibble()

# + magic_args="-o all_results" language="R"
# nrow(all_results)
# -

all_results

all_results.to_csv(output_file)
