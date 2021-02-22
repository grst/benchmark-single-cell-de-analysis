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

# import scanpy as sc
import pandas as pd
import numpy as np
from glob import glob
from pathlib import Path
import statsmodels.stats.multitest as multitest
from matplotlib_venn import venn3, venn2
import matplotlib.pyplot as plt

result_paths = [Path(x) for x in glob("../../data/results/*.csv")]

results = {p.stem: pd.read_csv(p) for p in result_paths}

list(results.keys())

tmp_scvi = results["scVI"].loc[
    :, ["Unnamed: 0", "proba_not_de", "lfc_mean", "comparison"]
]
tmp_scvi.columns = ["gene_symbol", "pvalue", "lfc", "comparison"]
results["scVI"] = tmp_scvi

tmp_diffxpy = results["diffxpy0.6.13"].loc[:, ["gene", "pval", "log2fc", "comparison"]]
tmp_diffxpy["comparison"] = tmp_diffxpy["comparison"].str.replace(".", "", regex=False)
tmp_diffxpy.columns = ["gene_symbol", "pvalue", "lfc", "comparison"]
results["diffxpy0.6.13"] = tmp_diffxpy

tmp_diffxpy7 = results["diffxpy0.7.4"].loc[:, ["gene", "pval", "log2fc", "comparison"]]
tmp_diffxpy7["comparison"] = tmp_diffxpy7["comparison"].str.replace(
    ".", "", regex=False
)
tmp_diffxpy7.columns = ["gene_symbol", "pvalue", "lfc", "comparison"]
results["diffxpy0.7.4"] = tmp_diffxpy7

tmp_sm = results["statsmodels_glm"].loc[
    :, ["gene", "pvalue", "log2_fold_change", "comparison"]
]
tmp_sm["comparison"] = tmp_sm["comparison"].str.replace(".", "", regex=False)
tmp_sm.loc[np.isnan(tmp_sm["pvalue"]), "pvalue"] = 1
tmp_sm.columns = ["gene_symbol", "pvalue", "lfc", "comparison"]
results["statsmodels_glm"] = tmp_sm

tmp_edger = results["edger"].loc[:, ["gene_symbol", "PValue", "logFC", "cluster"]]
tmp_edger["cluster"] = [
    f"{x.replace('leiden', '').strip()} vs Rest" for x in tmp_edger["cluster"]
]
tmp_edger.columns = ["gene_symbol", "pvalue", "lfc", "comparison"]
results["edger"] = tmp_edger

for key, tmp_results in results.items():
    if key == "scVI":
        tmp_results["fdr"] = tmp_results["pvalue"]
    else:
        tmp_results["fdr"] = multitest.fdrcorrection(
            tmp_results["pvalue"].values, method="indep"
        )[1]

for k, tmp_results in results.items():
    tmp_results.set_index("gene_symbol", inplace=True)


# ### Compare (top $n$)

def venn(keys, tmp_dict):
    if len(keys) == 3:
        return venn3([tmp_dict[k] for k in keys], set_labels=keys)
    elif len(keys) == 2:
        return venn2([tmp_dict[k] for k in keys], set_labels=keys)
    else:
        raise ValueError("only for 2 or 3 sets")


comparisions = list(results["edger"]["comparison"].unique())

comparisions

comparison = "12 vs Rest"

results_comparision = {
    k: df.loc[lambda x: x["comparison"] == comparison, :] for k, df in results.items()
}

top_genes = {
    method: set(
        tmp_results.sort_values("pvalue")
        .loc[lambda x: x["comparison"] == comparison, :]
        .loc[lambda x: x["lfc"] >= 1, :]
        .index.values[:100]
    )
    for method, tmp_results in results.items()
}

venn(["edger", "diffxpy0.7.4"], top_genes)

venn(["edger", "diffxpy0.6.13"], top_genes)

venn(["edger", "statsmodels_glm"], top_genes)

venn(["edger", "scVI"], top_genes)

venn(["scVI", "statsmodels_glm"], top_genes)

# ### Compare (all siginificant)

sig_genes = {
    method: set(
        tmp_results.sort_values("pvalue")
        .loc[lambda x: x["comparison"] == comparison, :]
        .loc[lambda x: x["fdr"] < 0.05, :]
        .index
    )
    for method, tmp_results in results.items()
}
for k, v in sig_genes.items():
    print(f"{k}: {len(v)}")

results["scVI"]

venn(["edger", "diffxpy0.7.4"], sig_genes)

venn(["edger", "diffxpy0.6.13"], sig_genes)

venn(["edger", "scVI"], sig_genes)

venn(["edger", "statsmodels_glm"], sig_genes)

venn(["scVI", "statsmodels_glm"], sig_genes)

venn(["diffxpy0.7.4", "statsmodels_glm"], sig_genes)

genes = list(results["edger"].index.values)

plt.scatter(
    results_comparision["edger"].loc[genes, "lfc"],
    results_comparision["diffxpy0.7.4"].loc[genes, "lfc"],
    s=1,
    alpha=0.05,
)

plt.scatter(
    results_comparision["diffxpy0.7.4"].loc[genes, "lfc"],
    results_comparision["diffxpy0.6.13"].loc[genes, "lfc"],
    s=1,
    alpha=0.05,
)

plt.scatter(
    results_comparision["statsmodels_glm"].loc[genes, "lfc"],
    results_comparision["diffxpy0.6.13"].loc[genes, "lfc"],
    s=1,
    alpha=0.05,
)

plt.scatter(
    results_comparision["edger"].loc[genes, "lfc"],
    results_comparision["statsmodels_glm"].loc[genes, "lfc"],
    s=1,
    alpha=0.05,
)

plt.scatter(
    results_comparision["edger"].loc[genes, "lfc"],
    results_comparision["scVI"].loc[genes, "lfc"],
    s=1,
    alpha=0.05,
)

# ### scatter pvalues

plt.scatter(
    -np.log10(results_comparision["edger"].loc[genes, "pvalue"]),
    -np.log10(results_comparision["diffxpy0.7.4"].loc[genes, "pvalue"]),
    s=1,
    alpha=0.05,
)

plt.scatter(
    -np.log10(results_comparision["edger"].loc[genes, "pvalue"]),
    -np.log10(results_comparision["statsmodels_glm"].loc[genes, "pvalue"]),
    s=1,
    alpha=0.05,
)

plt.scatter(
    -np.log10(results_comparision["diffxpy0.6.13"].loc[genes, "pvalue"]),
    -np.log10(results_comparision["statsmodels_glm"].loc[genes, "pvalue"]),
    s=1,
    alpha=0.05,
)

plt.scatter(
    -np.log10(results_comparision["edger"].loc[genes, "pvalue"]),
    -np.log10(results_comparision["scVI"].loc[genes, "pvalue"]),
    s=1,
    alpha=0.05,
)






