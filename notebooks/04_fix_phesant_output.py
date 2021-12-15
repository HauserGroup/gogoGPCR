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

# %%
import pandas as pd
import subprocess

# %%
PHENOTYPE="metabolic"
TRAIT="BT"

# %%
if TRAIT == "BT":
    file = f"/mnt/project/Data/phenotypes/{PHENOTYPE}.{TRAIT}.phesant/{PHENOTYPE}.phesantdata-binary-all.txt"
    
    df = pd.read_csv(file, sep=",", index_col="userID")
    
    df.columns = [s.replace("#", "_") for s in df.columns]
    
    for col in df.columns:
        
        assert df.loc[:, col].nunique() == 2, "Uh oh, column is not binary"
        
        df.loc[:, col] = df.loc[:, col].astype("Int64")
        
        print(col, "\n",df.loc[:, col].value_counts(), "\n", sep = "")

if TRAIT == "QT":
    file = f"/mnt/project/Data/phenotypes/{PHENOTYPE}.{TRAIT}.phesant/{PHENOTYPE}.phesantdata-cont-all.txt"
    
    df = pd.read_csv(file, sep=",", index_col="userID")
    
    df.columns = [s.replace("#", "_") for s in df.columns]
    
    for col in df.columns:
        print(df.loc[:, col].describe())

df.insert(0, "FID", df.index)
df.insert(1, "IID", df.index)

# %%
samp = pd.read_csv("../tmp/samples_in_200k_exomes.csv")
df = df.loc[samp.column,:]
df.head()

# %%
df.to_csv(
    f"../tmp/{PHENOTYPE}.{TRAIT}.final.tsv",
    sep="\t",
    na_rep="NA",
    index=False,
)

# %%
subprocess.run(
    ["dx", "upload", f"../tmp/{PHENOTYPE}.{TRAIT}.final.tsv", "--path", "Data/phenotypes/"],
    check=True,
    shell=False,
)
