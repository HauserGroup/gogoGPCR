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
import dxpy
import hail as hl

import pyspark
import tomli
import subprocess
from matrixtables import *
from utils import get_stats
from datetime import datetime
from pprint import pprint

Path("../tmp").resolve().mkdir(parents=True, exist_ok=True)

# %%
# Parameters
with open("../config.toml", "rb") as f:
    conf = tomli.load(f)

IMPORT = conf["IMPORT"]
NAME = conf["NAME"]
VCF_VERSION = IMPORT["VCF_VERSION"]
REFERENCE_GENOME = conf["REFERENCE_GENOME"]
DATABASE = IMPORT["DATABASE"]

LOG_FILE = (
    Path(IMPORT["LOG_DIR"], f"{NAME}_{datetime.now().strftime('%H%M')}.log")
    .resolve()
    .__str__()
)

GENE_FILE = Path(IMPORT["GENE_FILE"]).resolve().__str__()

with open(GENE_FILE, "r") as file:
    genes = file.read().splitlines()

if NAME == "NONE":
    NAME = genes[0]

# %%
# Spark and Hail

sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

mt_database = dxpy.find_one_data_object(name=DATABASE)["id"]
hl.init(sc=sc, default_reference=REFERENCE_GENOME, log=LOG_FILE)

# %%
STAGE = "QC1"
READ_PATH = "dnax://" + mt_database + f"/{NAME}.{STAGE}.mt"

mt = hl.read_matrix_table(READ_PATH)

v, s = mt.count()
pprint(f"{v} variants and {s} samples after reading matrixtable")

# %%
# Bespoke annotation table of interesting variants. In this case, it contains variants
# and Category columns

ht = hl.import_table(
    "file:" + "/mnt/project/Data/annotations/MC4R.tsv", impute=True, quote='"'
)
ht = ht.annotate(Variants=ht.Variants.strip())
ht = ht.key_by("Variants")

# %%
ht.show()

# %%
mt = mt.annotate_rows(labels=ht[mt.protCons].Category)

# %%
# Generate stats and summary of variants of interest

stats, intr = get_stats(mt)
stats.show(-1)

# %%
# Save QC info for variants of interest

intr.export(f"/tmp/{NAME}_QC1.tsv")

subprocess.run(
    ["hadoop", "fs", "-get", f"/tmp/{NAME}_QC1.tsv", f"../tmp/{NAME}_QC1.tsv"],
    check=True,
    shell=False,
)

# %%
# Checkpoint MT to database for loading into QC2 and further QC
STAGE = "LABELLED"
WRITE_PATH = "dnax://" + mt_database + f"/{NAME}.{STAGE}.mt"

mt.write(WRITE_PATH, overwrite=True)
