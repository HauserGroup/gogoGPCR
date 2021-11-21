# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Create hard sample filters
#
# For analysis, we need a cohort of samples with minimal population structure, minimal relatedness and without a few rare sources of error. This notebook generates a list of samples to remove from analysis in order to create such a cohort. We start out by importing stuff, initialising pyspark, setting various parameters from the configuration file, initialising Hail, and loading the participant dataset.

# %%
from distutils.version import LooseVersion
from pathlib import Path
import subprocess

import dxdata
import dxpy
import hail as hl
import pyspark
import tomli

from utils import fields_for_id

# %%
# Initialise Spark

sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

# %% tags=["parameters"]
# Set configurations

with open("../config.toml", "rb") as f:
    conf = tomli.load(f)

RAW_REL_FILE = conf["SAMPLE_QC"]["UKB_REL_DAT_FILE"]
FINAL_FILTER_FILE = conf["SAMPLE_QC"]["SAMPLE_FILTER_FILE"]

MAX_KINSHIP = conf["SAMPLE_QC"]["MAX_KINSHIP"]

LOG_FILE = Path(conf["IMPORT"]["LOG_DIR"], f"sample_filters.log").resolve().__str__()
TMP_DIR = Path(conf["EXPORT"]["TMP_DIR"])
DATA_DIR = Path(conf["SAMPLE_QC"]["DATA_DIR"])


# %%
# Initialise Hail

hl.init(sc=sc, default_reference="GRCh38", log=LOG_FILE)

# %%
# Load participant dataset

dispensed_database_name = dxpy.find_one_data_object(
    classname="database", name="app*", folder="/", name_mode="glob", describe=True
)["describe"]["name"]
dispensed_dataset_id = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob"
)["id"]

dataset = dxdata.load_dataset(id=dispensed_dataset_id)
participant = dataset["participant"]

# %% [markdown]
# # Filtering non-Caucasians and rare errors
# We filter out non-Caucasians (22006), outliers for heterozygosity or missing rate (22027), sex chromosome aneuploidy (22019) or genetic kinship to other participants (22021, UKB defined).

# %%
# Find relevant field names

fields = ["22027", "22019", "22006", "22021"]
field_names = [
    fields_for_id(id) for id in fields
]  # fields_for_id("22027") + fields_for_id("22019") + fields_for_id("22006") + fields_for_id("22021")
field_names = ["eid"] + [field.name for fields in field_names for field in fields]

# %%
# Retrieve dataframe

df = participant.retrieve_fields(
    names=field_names, engine=dxdata.connect(), coding_values="replace"
)

# %%
df.show(5, truncate=False)

# %%
# Use hard filters

df = df.filter(
    df.p22006.isNull()
    | (~df.p22027.isNull())
    | (~df.p22019.isNull())
    | (df.p22021 == "Participant excluded from kinship inference process")
    | (df.p22021 == "Ten or more third-degree relatives identified")
)
filtered_samples_to_remove = hl.Table.from_spark(df.select("eid")).key_by("eid")
print(f"Samples to be filtered: {filtered_samples_to_remove.count()}")

# %% [markdown]
# # Filter related samples
# UK Biobank provides a list of genetically related individuals (KING) called 'ukb_rel.dat' which contains a kinship coefficient between pairs of individuals. Here, we remove any sample with a closer than 3rd degree relative (kinship > 0.088) and which is not already filtered out in the previous step. We then use Hail to create a maximal independent set of individuals by removing the smallest amount of related individuals. This is finally combined with the previously filtered samples to give the final list of samples to remove from the analysis.

# %%
# Import related table, remove any individual already sampled and keep those with kinship > 0.088

rel = hl.import_table(
    "file:" + "/mnt/project/" + RAW_REL_FILE,
    delimiter=" ",
    impute=True,
    types={"ID1": "str", "ID2": "str"},
)

rel = (
    rel.key_by("ID2")
    .anti_join(filtered_samples_to_remove)
    .key_by("ID1")
    .anti_join(filtered_samples_to_remove)
)

rel = rel.filter(rel.Kinship > MAX_KINSHIP, keep=True)

print(
    f"Related samples not already in filter and low kinship coefficient: {rel.count()}"
)

# %%
# Find maximal independent set

related_samples_to_remove = (
    hl.maximal_independent_set(
        i=rel.ID1,
        j=rel.ID2,
        keep=False,
    )
    .rename({"node": "eid"})
    .key_by("eid")
)

print(
    f"Samples to remove to create independent set: {related_samples_to_remove.count()}"
)

# %%
# Join the two sets of samples to remove

final = related_samples_to_remove.join(filtered_samples_to_remove, how="outer")
print(f"Final number of samples to remove: {final.count()}")

# %%
# Export list

FILTER_PATH = (TMP_DIR / FINAL_FILTER_FILE).resolve().__str__()
PROCESSED_DIR = (DATA_DIR.parents[0].stem / Path(DATA_DIR.stem)).__str__() + "/"

final.export("file:" + FILTER_PATH)

# %%
# Upload to project

subprocess.run(
    ["dx", "upload", FILTER_PATH, "--path", PROCESSED_DIR], check=True, shell=False
)
