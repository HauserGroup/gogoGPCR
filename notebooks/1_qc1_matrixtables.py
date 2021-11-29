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
FILTER_FILE = (
    Path(conf["SAMPLE_QC"]["DATA_DIR"], conf["SAMPLE_QC"]["SAMPLE_FILTER_FILE"])
    .resolve()
    .__str__()
)

VCF_DIR = Path(IMPORT["VCF_DIR"]).resolve().__str__()

DOWNSAMPLE_P = IMPORT.get("DOWNSAMPLE_P", None)

SNV_ONLY = conf["ANNOTATE"]["SNV_ONLY"]
USE_VEP = conf["ANNOTATE"]["USE_VEP"]
MISSENSE_ONLY = conf["ANNOTATE"]["MISSENSE_ONLY"]

VEP_JSON = Path(conf["ANNOTATE"]["VEP_JSON"]).resolve().__str__()

ANNOTATION_DIR = conf["ANNOTATE"]["ANNOTATION_DIR"]

MIN_DP = conf["ENTRY_QC"]["MIN_DP"]
MIN_GQ = conf["ENTRY_QC"]["MIN_GQ"]
MIN_PL = conf["ENTRY_QC"]["MIN_PL"]

MIN_P_HWE = conf["VARIANT_QC"]["MIN_P_HWE"]
MIN_VAR_GQ = conf["VARIANT_QC"]["MIN_VAR_GQ"]

MIN_CALL_RATE = conf["SAMPLE_QC"]["MIN_CALL_RATE"]
MIN_MEAN_DP = conf["SAMPLE_QC"]["MIN_MEAN_DP"]
MIN_MEAN_GQ = conf["SAMPLE_QC"]["MIN_MEAN_GQ"]

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
with open(GENE_FILE, "r") as file:
    genes = file.read().splitlines()
    
if NAME == "NONE":
    NAME = genes[0]

mapping = pd.read_csv(MAP_FILE, sep="\t").set_index("HGNC", drop=False)

# %%
# Import
mt = import_mt(genes, mapping, vcf_dir=VCF_DIR, vcf_version=VCF_VERSION).key_rows_by(
    "locus", "alleles"
)  # .checkpoint(checkpoint_file)

v, s = mt.count()
pprint(f"{v} variants and {s} samples after import")

# %%
# Checkpoint
stage = "raw"
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

mt.write(WRITE_PATH, overwrite = True)

# %%
ANNOTATION_FILE = Path(ANNOTATION_DIR, f"{NAME}.tsv").resolve().__str__()
mt = annotate_mt(mt=mt, gene=NAME, annotations=ANNOTATION_FILE)

interesting = mt.filter_rows(
    (hl.is_defined(mt.annotations)) & (hl.agg.any(mt.GT.is_non_ref()))
).count_rows()
pprint(f"{interesting} annotated variants found before QC")

# %%
# Checkpoint
stage = "QC1"
checkpoint_file = f"/tmp/{NAME}.{stage}.cp.mt"

mt = mt.checkpoint(checkpoint_file, overwrite=True)
# show_stats(mt)

# %%
# Withdrawn
mt = mt.filter_cols(~mt.s.startswith("W"))

print(f"Samples remaining after removing withdrawn participants: {mt.count_cols()} ")

# %%
# Filter samples
samples_to_remove = hl.import_table("file:" + FILTER_FILE, key="eid")
mt = mt.anti_join_cols(samples_to_remove)
print(f"Samples remaining after removing related samples: {mt.count_cols()} ")

# %%
# Sample QC
mt = sample_QC_mt(mt, MIN_CALL_RATE, MIN_MEAN_DP, MIN_MEAN_GQ)

print(f"Samples remaining after QC: {mt.count_cols()} ")

# %%
# Variant QC
mt = variant_QC_mt(mt, MIN_P_HWE, MIN_VAR_GQ)

interesting = mt.filter_rows(
    (hl.is_defined(mt.annotations)) & (hl.agg.any(mt.GT.is_non_ref()))
).count_rows()
print(
    f"{mt.count_rows()} variants remaining after QC of which {interesting} are annotated"
)

# %%
# Genotype GQ
mt = genotype_filter_mt(mt, MIN_DP, MIN_GQ, True)

missing = mt.aggregate_entries(hl.agg.sum(~hl.is_defined(mt.GT)))
pprint(f"{missing} missing or filtered entries after Call QC")

# %%
# Checkpoint
stage = "QC2"
checkpoint_file = f"/tmp/{GENE}.{stage}.cp.mt"

mt = mt.checkpoint(checkpoint_file, overwrite=True)
show_stats(mt)

# %%
# BGEN
write_bgen(mt, "file:" + BGEN_FILE)

# %%
# ANNOTATIONS

mt = add_varid(mt)

annotations = (
    mt.select_rows(
        varid=mt.varid,
        gene=mt.vep.transcript_consequences.gene_symbol[0],
        annotation=mt.annotation,
    )
    .rows()
    .key_by("varid")
    .drop("locus")
    .drop("alleles")
)
annotations.export("file:" + ANNOTATIONS_FILE, header=False)

# %%
# SETLIST
position = mt.aggregate_rows(hl.agg.min(mt.locus.position))
names = mt.varid.collect()
names_str = ",".join(names)

line = f"{mt.vep.transcript_consequences.gene_symbol[0].collect()[0]}\t{mt.locus.contig.collect()[0]}\t{position}\t{names_str}"

with open(SETLIST_FILE, "w") as f:
    f.write(line)

# %%
bgen_file = BGEN_FILE + ".bgen"
sample_file = BGEN_FILE + ".sample"

# subprocess.run(["dx", "upload", bgen_file, sample_file, ANNOTATIONS_FILE, SETLIST_FILE, "--path", "/data/burden/"], check = True, shell = False)

# %%
sample = mt.select_cols(ID_1=mt.s, ID_2=mt.s, missing=0)

# %%
sample.cols().show()

# %%
# STAGE = "final"
# WRITE_PATH = "dnax://" + mt_database + f"/{GENE}.{STAGE}.mt"

# mt.write(WRITE_PATH, overwrite = True)
show_stats(mt)

# STAGE = "final"
# WRITE_PATH = "dnax://" + mt_database + f"/{GENE}.{STAGE}.mt"

# mt = hl.read_matrix_table(WRITE_PATH)
