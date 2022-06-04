# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.6
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Forest plots for quantitative and binary regenie results
# Forest plots can be generated for both binary and quantitative results.
# See data/examples for formatted regenie output files
# Everything here is hacky and use at your own peril

# %%
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from results import pheno_search, pval_stars, plot_BT, plot_QT

Path("../tmp").resolve().mkdir(parents=True, exist_ok=True)

# %%
# Flags

GENE = "MC4R"
TRAIT = "BT"
PHENOTYPE = "metabolic"
AAF = None
# %%
# Results files

files = [
    f"file:/mnt/project/Data/results/{PHENOTYPE}.{TRAIT}.{GENE}/{file}"
    for file in os.listdir(f"/mnt/project/Data/results/{PHENOTYPE}.{TRAIT}.{GENE}")
    if file.endswith(".regenie")
]

# Field codings
ukb_coding = pd.read_csv(
    "/opt/notebooks/gogoGPCR/data/misc/Data_Dictionary_Showcase.csv",
    error_bad_lines=False,
    warn_bad_lines=False,
    quotechar='"',
    usecols=["FieldID", "Field"],
)

custom_coding = pd.read_csv(
    "/opt/notebooks/gogoGPCR/data/misc/Data_Dictionary_Custom.csv",
    error_bad_lines=False,
    warn_bad_lines=False,
    quotechar='"',
)

# %%
# Load raw DF

df_raw = pd.read_csv(files[0], delimiter=" ", header="infer", comment="#").assign(
    SOURCE=os.path.basename(files[0])
)
df_raw = pd.concat(
    [df_raw]
    + [
        pd.read_csv(fp, delimiter=" ", comment="#").assign(SOURCE=os.path.basename(fp))
        for fp in files[1:]
    ],
    axis=0,
)

# %%
# Fix common fields
df = df_raw
df.loc[:, "GENE"] = df.ID.apply(lambda x: x.split(".")[0])
df.loc[:, "MASK"] = df.ALLELE1.apply(lambda x: x.split(".", maxsplit=2)[0])
df.loc[:, "AAF"] = df.ALLELE1.apply(lambda x: x.split(".", maxsplit=1)[-1])
df.loc[:, "TRAIT"] = TRAIT  # df.FILE.apply(lambda x: x[0])
df.loc[:, "PHENO"] = df.SOURCE.apply(lambda x: x.split("_")[-1].split(".")[0])
df = df.drop(["ID", "ALLELE0", "ALLELE1", "EXTRA", "SOURCE", "TEST"], axis=1)

# Filter AAF
if not AAF:
    AAF = df.AAF.max()

df = df.loc[df.AAF == AAF, :]

# Sanity check
df.head()

# %%
# Filters
bt = df.TRAIT == "BT"
qt = df.TRAIT == "QT"

# Fix Binary Traits
df.loc[bt, "OR"] = np.exp(df.loc[bt, "BETA"])
df.loc[bt, "OR_up"] = np.exp(df.loc[bt, "BETA"] + df.loc[bt, "SE"])
df.loc[bt, "OR_low"] = np.exp(df.loc[bt, "BETA"] - df.loc[bt, "SE"])
df.loc[bt, "OR_up_lim"] = df.loc[bt, "OR_up"] - df.loc[bt, "OR"]
df.loc[bt, "OR_low_lim"] = df.loc[bt, "OR"] - df.loc[bt, "OR_low"]

# Fix Quantitative Traits
df.loc[qt, "BETA_up_lim"] = df.loc[qt, "BETA"] + df.loc[qt, "SE"]
df.loc[qt, "BETA_low_lim"] = df.loc[qt, "BETA"] - df.loc[qt, "SE"]

# Final fixes
df.loc[:, "Phenotype"] = df.PHENO.apply(
    lambda x: pheno_search(x, ukb_coding, custom_coding).replace('"', "").strip()
)
df.loc[:, "pval"] = np.power(10, -df["LOG10P"])
df.loc[:, "pval_stars"] = df["pval"].apply(lambda x: pval_stars(x))
df.loc[:, "N_pos"] = (2 * df["N"] * df["A1FREQ"]).astype(int)

# Singletons
df = df.loc[df.AAF != "singleton", :]
print(len(df))
df.head()

# %%
phenos_to_remove = []

plt_df = (
    df.loc[(df.GENE == GENE)]
    .sort_values(by=["Phenotype", "AAF"], ascending=[True, False])  # , "AAF"
    .groupby(["Phenotype", "MASK"])
    .first()
    .reset_index()
)

plt_df = plt_df.loc[~plt_df.Phenotype.astype(str).isin(phenos_to_remove), :]

effect = {"BT": "OR", "QT": "BETA"}[TRAIT]

group_by_mean = (
    pd.DataFrame({"mean": plt_df.groupby(["Phenotype"]).agg("mean")[effect]})
    .sort_values(by="mean", ascending=False)
    .reset_index()
)

sorter = group_by_mean.Phenotype.tolist()

plt_df.loc[:, "Phenotype"] = plt_df.loc[:, "Phenotype"].astype("category")
plt_df.loc[:, "Phenotypes"] = plt_df.loc[:, "Phenotype"].cat.set_categories(sorter)

plt_df = plt_df.sort_values(
    by=["Phenotype", "MASK"], ascending=[True, False]
).reset_index(drop=True)

phenotypes = plt_df.Phenotype.unique()

print(len(plt_df))
plt_df.head()


# %%
plot = plot_BT(
    plt_df, title=f"{GENE} Binary {PHENOTYPE.capitalize()} Traits", xlim=[0, 10]
)

plt.savefig(
    f"/opt/notebooks/gogoGPCR/tmp/{GENE}.{TRAIT}.svg",
    dpi=600,
    bbox_inches="tight",
    format="svg",
)


# %%

plot = plot_QT(
    plt_df,
    title=f"{GENE} Quantitative {PHENOTYPE.capitalize()} Traits",
    xlim=[-1, 1],
    height=12,
)


plt.savefig(
    f"/opt/notebooks/gogoGPCR/tmp/{GENE}.{TRAIT}.svg",
    dpi=600,
    bbox_inches="tight",
    format="svg",
)

# %%
plt_df.to_csv(f"/opt/notebooks/gogoGPCR/tmp/{GENE}.{TRAIT}.csv")

# %%
