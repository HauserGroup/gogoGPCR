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
# # Create a MatrixTable and QC the hell out of it
# ## Import stuff and set your parameters
# First, we import necessary libraries and configurations from config.toml. Then we initialise Spark and Hail.

# %%
import subprocess
from datetime import datetime
from distutils.version import LooseVersion
from functools import partial
from pathlib import Path
from pprint import pprint

import dxdata
import dxpy
import hail as hl
import pandas as pd
import pyspark
import tomli
from matrixtables import *
from utils import show_stats

Path("../tmp").resolve().mkdir(parents=True, exist_ok=True)

# %% tags=["parameters"]
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

MAP_FILE = Path(IMPORT["MAPPING_FILE"]).resolve().__str__()
INT_FILE = Path(IMPORT["INTERVAL_FILE"]).resolve().__str__()
GENE_FILE = Path(IMPORT["GENE_FILE"]).resolve().__str__()

with open(GENE_FILE, "r") as file:
    GENES = file.read().splitlines()

    if NAME == "NONE":
        NAME = GENES[0]


VCF_DIR = Path(IMPORT["VCF_DIR"]).resolve().__str__()

DOWNSAMPLE_P = IMPORT.get("DOWNSAMPLE_P", None)

SNV_ONLY = conf["ANNOTATE"]["SNV_ONLY"]
USE_VEP = conf["ANNOTATE"]["USE_VEP"]
MISSENSE_ONLY = conf["ANNOTATE"]["MISSENSE_ONLY"]

VEP_JSON = Path(conf["ANNOTATE"]["VEP_JSON"]).resolve().__str__()

ANNOTATION_DIR = conf["ANNOTATE"]["ANNOTATION_DIR"]


TMP_DIR = conf["EXPORT"]["TMP_DIR"]

BGEN_FILE = Path(TMP_DIR, f"{NAME}").resolve().__str__()
ANNOTATIONS_FILE = Path(TMP_DIR, f"{NAME}.annotations").resolve().__str__()
SETLIST_FILE = Path(TMP_DIR, f"{NAME}.setlist").resolve().__str__()


# %%
# Spark and Hail

sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

try:
    mt_database = dxpy.find_one_data_object(name=DATABASE)["id"]
except Exception as e:
    spark.sql(f"CREATE DATABASE {DATABASE} LOCATION  'dnax://'")
    mt_database = dxpy.find_one_data_object(name=DATABASE)["id"]

# this breaks export_bgen for now
# hl.init(sc=sc, default_reference=REFERENCE_GENOME, log=LOG_FILE, tmp_dir=f'dnax://{mt_database}/tmp/')

hl.init(sc=sc, default_reference=REFERENCE_GENOME, log=LOG_FILE)

# %%
# Import
mapping = pd.read_csv(MAP_FILE, sep="\t").set_index("HGNC", drop=False)

mt = import_mt(GENES, mapping, vcf_dir=VCF_DIR, vcf_version=VCF_VERSION).key_rows_by(
    "locus", "alleles"
)  # .checkpoint(checkpoint_file)

v, s = mt.count()
pprint(f"{v} variants and {s} samples after import")

# %%
# Checkpoint
stage = "RAW"
checkpoint_file = f"/tmp/{NAME}.{stage}.cp.mt"

mt = mt.checkpoint(checkpoint_file, overwrite=True)

# %%
# Downsample
if DOWNSAMPLE_P is not None:
    mt = downsample_mt(mt, DOWNSAMPLE_P)

    pprint(f"{mt.count_cols()} samples after downsampling")

# %%
# Interval QC
mt = interval_qc_mt(mt, "file:" + INT_FILE)

pprint(f"{mt.count_rows()} variants after interval filtering")

# %%
# Split multi
mt = mt.filter_rows(mt.alleles.length() <= 6)
mt = smart_split_multi_mt(mt)

pprint(f"{mt.count_rows()} variants with not more than 6 alleles after splitting")

# %%
if USE_VEP:
    mt = hl.vep(mt, "file:" + VEP_JSON)

    is_MANE = mt.aggregate_rows(
        hl.agg.all(hl.is_defined(mt.vep.transcript_consequences.mane_select))
    )
    assert is_MANE, "Selected transcript may not be MANE Select. Check manually."

    mt = mt.annotate_rows(
        protCons=mt.vep.transcript_consequences.amino_acids[0].split("/")[0]
        + hl.str(mt.vep.transcript_consequences.protein_end[0])
        + mt.vep.transcript_consequences.amino_acids[0].split("/")[-1]
    )


# %%
STAGE = "QC1"
WRITE_PATH = "dnax://" + mt_database + f"/{NAME}.{STAGE}.mt"

mt.write(WRITE_PATH, overwrite=True)
