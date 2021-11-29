# ---
# jupyter:
#   jupytext:
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

# %%
import dxdata
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
ht = hl.import_table("file:" + "/mnt/project/data/annotations/DRD2.tsv", impute=True).key_by("AA consequence")

# %%
mt = mt.annotate_rows(annotations=ht[mt.protCons])

# %%
mt = mt.annotate_rows(
        Gi1=mt.annotations.number_of_impairments_Gi1 > 0,
        GoA=mt.annotations.number_of_impairments_GoA > 0,
        Gz=mt.annotations.number_of_impairments_Gz > 0,
    )

mt = mt.annotate_rows(
    labels=hl.case()
    .when(~mt.Gi1 & ~mt.GoA & ~mt.Gz, "WT")
    .when(mt.Gi1 & ~mt.GoA & ~mt.Gz, "Gi1")
    .when(~mt.Gi1 & mt.GoA & ~mt.Gz, "GoA")
    .when(~mt.Gi1 & ~mt.GoA & mt.Gz, "Gz")
    .when(mt.Gi1 & mt.GoA & ~mt.Gz, "Gi1_GoA")
    .when(mt.Gi1 & ~mt.GoA & mt.Gz, "Gi1_Gz")
    .when(~mt.Gi1 & mt.GoA & mt.Gz, "GoA_Gz")
    .when(mt.Gi1 & mt.GoA & mt.Gz, "Gi1_GoA_Gz")
    .or_missing()
)

# %%
stats, intr = get_stats(mt)
stats.show(-1)

# %%
intr.export(f"/tmp/{NAME}_QC1.tsv")

subprocess.run(
    ["hadoop", "fs", "-get", f"/tmp/{NAME}_QC1.tsv", f"../tmp/{NAME}_QC1.tsv"],
    check=True,
    shell=False,
)

# %%
STAGE = "LABELLED"
WRITE_PATH = "dnax://" + mt_database + f"/{NAME}.{STAGE}.mt"

mt.write(WRITE_PATH, overwrite = True)
