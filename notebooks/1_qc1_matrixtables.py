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
# # Initial quality control (QC1)
# We import the WES VCF files as MatrixTables and perform QC common to all workflows

# %%
# Imports

import subprocess
from datetime import datetime
from distutils.version import LooseVersion
from functools import partial
from pathlib import Path
from pprint import pprint

import dxpy
import hail as hl
import pandas as pd
import pyspark
import tomli
from matrixtables import *
from utils import get_stats

Path("../tmp").resolve().mkdir(parents=True, exist_ok=True)

# %% tags=["parameters"]
# Parameters

# Config file
# As far as possible this contains all configurations
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
hl.init(sc=sc, default_reference=REFERENCE_GENOME, log=LOG_FILE, tmp_dir=f'dnax://{mt_database}/tmp/')

hl.init(sc=sc, default_reference=REFERENCE_GENOME, log=LOG_FILE)


# %% [markdown]
# ## Import a list of VCF files corresponding to genes from genes.txt.
# Note that hl.import_gvcfs needs a list of regions to import, here we use those from the provided genes
# See MAP_FILE for an example of necessary columns for import_mt. The one provided in here, in ../data,
# only contains GPCR genes
# %%
# Import MatrixTable
mapping = pd.read_csv(MAP_FILE, sep="\t").set_index("HGNC", drop=False)

mt = import_mt(GENES, mapping, vcf_dir=VCF_DIR, vcf_version=VCF_VERSION).key_rows_by(
    "locus", "alleles"
)

v, s = mt.count()
pprint(f"{v} variants and {s} samples after import")

# %%
# Checkpoint
# Hail likes checkpoints after major operations

stage = "RAW"
checkpoint_file = f"/tmp/{NAME}.{stage}.cp.mt"

mt = mt.checkpoint(checkpoint_file, overwrite=True)

# %%
NAME = "MC4R"
stage = "RAW"
mt = hl.read_matrix_table(f"/tmp/{NAME}.{stage}.cp.mt")

# %%
# Downsample
# If provided, samples can be downsampled for test purposes

if DOWNSAMPLE_P is not None:
    mt = downsample_mt(mt, DOWNSAMPLE_P)

    pprint(f"{mt.count_cols()} samples after downsampling")

# %%
# Interval QC
# Filter to only WES capture regions

mt = interval_qc_mt(mt, "file:" + INT_FILE)

pprint(f"{mt.count_rows()} variants after interval filtering")

# %%
# Split multi
# Create MatrixTable of only biallelic sites. smart_split_multi_mt may be over-engineered
# but should be marginally faster for large datasets

mt = mt.filter_rows(mt.alleles.length() <= 6)
mt = smart_split_multi_mt(mt)

pprint(f"{mt.count_rows()} variants with not more than 6 alleles after splitting")

# %%
# Annotate variants with VEP
# As of 104, MANE_SELECT is the canonical transcript and this function may break or may not

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
# Save QC'ed MT to database for loading into separate annotation notebook (e.g. 2a or 2b)
# Can also be used straight in QC2

STAGE = "QC1"
WRITE_PATH = f"/tmp/{NAME}.{STAGE}.mt"

mt.write(WRITE_PATH, overwrite=True)

# %%
