# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.10.0
#   kernelspec:
#     display_name: Python [conda env:.conda-2021-sc-de-benchmark-diffxpy07]
#     language: python
#     name: conda-env-.conda-2021-sc-de-benchmark-diffxpy07-py
# ---

from tqdm.contrib.concurrent import process_map
import itertools
import warnings
import contextlib
import scanpy as sc
from tqdm import tqdm
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError
from nxfvars import nxf

adata = sc.read_h5ad(nxf.input("adata", "../../data/adata-myeloid.h5ad"))
output_file = nxf.input("output_file", "/dev/null")


def test_gene(i, ad):
    data = pd.DataFrame().assign(
        gene=ad.X[:, i].todense().A1,
        dataset=ad.obs["dataset"].values,
        n_genes=ad.obs["n_genes"].values,
        group=ad.obs["leiden"].values,
    )
    f = open("/dev/null", "w")
    with contextlib.redirect_stdout(f):
        with contextlib.redirect_stderr(f):
            try:
                model = sm.GLM.from_formula(
                    "gene ~ 0+ C(group) + C(dataset) + n_genes",
                    data=data,
                    family=sm.families.NegativeBinomial(),
                ).fit()
                return model
            except (ValueError, np.linalg.linalg.LinAlgError, PerfectSeparationError):
                pass


model = test_gene(0, adata)

# %%time
models = process_map(
    test_gene,
    range(adata.shape[1]),
    itertools.repeat(adata),
    chunksize=1000,
    max_workers=int(nxf.task("cpus", "16")),
    tqdm_class=tqdm,
)

model = models[0]

model.summary()

tests = dict()
groups = adata.obs["leiden"].unique()
for group in groups:
    contrast = " + ".join([f"C(group)[{g}]" for g in groups[groups != group]])
    tests[group] = f"C(group)[{group}] - ({contrast})/{ len(groups) - 1}"

records = []
for group, test in tests.items():
    for i, model in enumerate(tqdm(models)):
        test = model.wald_test("C(group)[0] - (C(group)[8] + C(group)[10] + C(group)[12])/3") if model is not None else None
        gene = adata.X[:, i].todense().A1
        fold_change = np.mean(gene[adata.obs["leiden"] == group]) / np.mean(gene[adata.obs["leiden"] != group])
        records.append({
            "gene": adata.var_names[i],
            "comparison": f"{group} vs. Rest",
            "fold_change": fold_change,
            "log2_fold_change": np.log2(fold_change),
            "pvalue": test.pvalue if test is not None else np.nan
        })

res = pd.DataFrame.from_records(
    records
)

res.sort_values("pvalue").head(20)

res.to_csv(output_file)
