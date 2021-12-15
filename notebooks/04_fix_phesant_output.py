# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.4
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Fix the PHESANT output for regenie
# Again, more wrangling. For ease of use and further analysis, dataframes must contain exlusively quantitative or binary columns—not a mix


# %%
import pandas as pd
import subprocess

# %%
# Set FLAGS

PHENOTYPE = "metabolic"
TRAIT = "BT"

# %%
# Load output from PHESANT and fix headers

if TRAIT == "BT":
    file = f"/mnt/project/Data/phenotypes/{PHENOTYPE}.{TRAIT}.phesant/{PHENOTYPE}.phesantdata-binary-all.txt"

    df = pd.read_csv(file, sep=",", index_col="userID")

    df.columns = [s.replace("#", "_") for s in df.columns]

    for col in df.columns:

        assert df.loc[:, col].nunique() == 2, "Uh oh, column is not binary"

        df.loc[:, col] = df.loc[:, col].astype("Int64")

        print(col, "\n", df.loc[:, col].value_counts(), "\n", sep="")

if TRAIT == "QT":
    file = f"/mnt/project/Data/phenotypes/{PHENOTYPE}.{TRAIT}.phesant/{PHENOTYPE}.phesantdata-cont-all.txt"

    df = pd.read_csv(file, sep=",", index_col="userID")

    df.columns = [s.replace("#", "_") for s in df.columns]

    for col in df.columns:
        print(df.loc[:, col].describe())

# regenie expects these two first columns

df.insert(0, "FID", df.index)
df.insert(1, "IID", df.index)

# %%
# Save with tabs and upload

df.to_csv(
    f"../tmp/{PHENOTYPE}.{TRAIT}.final.tsv",
    sep="\t",
    na_rep="NA",
    index=False,
)


subprocess.run(
    [
        "dx",
        "upload",
        f"../tmp/{PHENOTYPE}.{TRAIT}.final.tsv",
        "--path",
        "Data/phenotypes/",
    ],
    check=True,
    shell=False,
)
