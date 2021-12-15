# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# # Generate psychiatric phenotypes dataframe
# An example of how to pull phenotypes and wrangle them into a form that PHESANT likes
# %%
# Imports

import re
import subprocess

import dxdata
import dxpy
import pandas as pd
import pyspark

from pyspark.sql import functions as F
from pyspark.conf import SparkConf
from pyspark.sql.types import StringType

from pathlib import Path

from phenotypes import *

Path("../tmp").resolve().mkdir(parents=True, exist_ok=True)

# %%
# Access the spark dataset

conf = SparkConf()
conf.set("autoBroadcastJoinThreshold", -1)

sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

# As recommended by RAP documentation
dispensed_database_name = dxpy.find_one_data_object(
    classname="database", name="app*", folder="/", name_mode="glob", describe=True
)["describe"]["name"]
dispensed_dataset_id = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob"
)["id"]

spark.sql("USE " + dispensed_database_name)

dataset = dxdata.load_dataset(id=dispensed_dataset_id)
participant = dataset["participant"]

# %% [markdown]
# We create a list of every field name in the dataset that may be of interest to us including the age and sex fields necessary for PHESANT
# %%
# Get all first occurrence fields for psychiatric phenotypes

first_occurences = list(
    participant.find_fields(
        lambda f: bool(re.match("^Date [F]\d{2} first reported", f.title))
    )
)

first_occurence_fields = [f.name for f in first_occurences]

# Age and sex fields

age_sex_fields = get_age_sex(participant, fields=["31", "21022"])

# Miscellaneous psychiatric fields

psych_fields = get_pheno_fields(
    participant,
    fields=[
        "2090",
        "2100",
        "20126",
    ],
)

# Concatenate the fields

field_names = concatenate(["eid"], age_sex_fields, first_occurence_fields, psych_fields)

# Create dataframe

df = participant.retrieve_fields(names=field_names, engine=dxdata.connect())

# %%
# Only keeo fields with at least 500 cases, any less is probably going to get wonky

min_cases = 500

print(f"Number of columns: {len(df.columns)}")
to_keep = get_columns_to_keep(df, min_cases)

df = df.select(*to_keep)

print(f"Number of columns with at least {min_cases} cases: {len(df.columns)}")

# %%
# Fix the headers for PHESANT and only include 200k samples with exome data

df = fix_colnames(df)
df = filter_to_200k(df)

# %%
df = df.drop(
    "x130836_0_0", "x130838_0_0", "x130842_0_0"
)  # these do not converge during Firth correction

# %%
# Get files from HDFS to local

df.write.csv("/tmp/phenos.tsv", sep="\t", header=True, emptyValue="NA")

subprocess.run(
    ["hadoop", "fs", "-rm", "/tmp/phenos.tsv/_SUCCESS"], check=True, shell=False
)
subprocess.run(
    ["hadoop", "fs", "-get", "/tmp/phenos.tsv", "../tmp/phenos.tsv"],
    check=True,
    shell=False,
)

# %%
# Fix files
# !sed -e '3,${/^xeid/d' -e '}' ../tmp/phenos.tsv/part* > ../tmp/psychiatric.raw.tsv

# %%
# Upload to project

subprocess.run(
    ["dx", "upload", "../tmp/psychiatric.raw.tsv", "--path", "Data/phenotypes/"],
    check=True,
    shell=False,
)
