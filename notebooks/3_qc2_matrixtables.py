# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

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
from utils import get_stats

Path("../tmp").resolve().mkdir(parents=True, exist_ok=True)

# %%
with open("../config.toml", "rb") as f:
    conf = tomli.load(f)

# BASICS
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

FILTER_FILE = (
    Path(conf["SAMPLE_QC"]["DATA_DIR"], conf["SAMPLE_QC"]["SAMPLE_FILTER_FILE"])
    .resolve()
    .__str__()
)

# GENES
GENE_FILE = Path(IMPORT["GENE_FILE"]).resolve().__str__()

with open(GENE_FILE, "r") as file:
    GENES = file.read().splitlines()
    
if NAME == "NONE":
    NAME = GENES[0]

# SAMPLE        
MIN_CALL_RATE = conf["SAMPLE_QC"]["MIN_CALL_RATE"]
MIN_MEAN_DP = conf["SAMPLE_QC"]["MIN_MEAN_DP"]
MIN_MEAN_GQ = conf["SAMPLE_QC"]["MIN_MEAN_GQ"]

# VARIANT
MIN_P_HWE = conf["VARIANT_QC"]["MIN_P_HWE"]
MIN_VAR_GQ = conf["VARIANT_QC"]["MIN_VAR_GQ"]

# GENOTYPE
MIN_DP = conf["ENTRY_QC"]["MIN_DP"]
MIN_GQ = conf["ENTRY_QC"]["MIN_GQ"]

# EXPORT
TMP_DIR = conf["EXPORT"]["TMP_DIR"]

BGEN_FILE = Path(TMP_DIR, f"{NAME}").resolve().__str__()
ANNOTATIONS_FILE = Path(TMP_DIR, f"{NAME}.annotations").resolve().__str__()
SETLIST_FILE = Path(TMP_DIR, f"{NAME}.setlist").resolve().__str__()


# %%
# Spark and Hail
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

mt_database = dxpy.find_one_data_object(name=DATABASE)["id"]
hl.init(sc=sc, default_reference=REFERENCE_GENOME, log=LOG_FILE)

# %%
STAGE = "LABELLED"
READ_PATH = "dnax://" + mt_database + f"/{NAME}.{STAGE}.mt"

mt = hl.read_matrix_table(READ_PATH)

interesting = mt.filter_rows(
    (hl.is_defined(mt.labels)) & (hl.agg.any(mt.GT.is_non_ref()))
).count_rows()

pprint(f"{interesting} annotated variants found before QC")

# %%
# Withdrawn
mt = mt.filter_cols(~mt.s.startswith("W"))

pprint(f"Samples remaining after removing withdrawn participants: {mt.count_cols()} ")

# %%
# Filter samples
samples_to_remove = hl.import_table("file:" + FILTER_FILE, key="eid")
mt = mt.anti_join_cols(samples_to_remove)

pprint(f"Samples remaining after hard filtering samples: {mt.count_cols()} ")

# %%
# Sample QC
# May need adjustment if too many samples are removed by default settings
#MIN_MEAN_DP = 15
#MIN_MEAN_GQ = 48.5

mt = sample_QC_mt(mt, MIN_CALL_RATE, MIN_MEAN_DP, MIN_MEAN_GQ)

pprint(f"Samples remaining after QC: {mt.count_cols()} ")

# %%
# Variant QC
mt = variant_QC_mt(mt, MIN_P_HWE, MIN_VAR_GQ)

interesting = mt.filter_rows(
    (hl.is_defined(mt.labels)) & (hl.agg.any(mt.GT.is_non_ref()))
).count_rows()

pprint(
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
checkpoint_file = f"/tmp/{NAME}.{stage}.cp.mt"

mt = mt.checkpoint(checkpoint_file, overwrite=True)

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
        annotation=mt.labels,
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

subprocess.run(["dx", "upload", bgen_file, sample_file, ANNOTATIONS_FILE, SETLIST_FILE, "--path", "/Data/burden/"], check = True, shell = False)

# %%
STAGE = "FINAL"
WRITE_PATH = "dnax://" + mt_database + f"/{NAME}.{STAGE}.mt"

mt.write(WRITE_PATH, overwrite = True)

# %%
stats, intr = get_stats(mt)
stats.show(-1)

# %%
intr.export(f"/tmp/{NAME}_QC2.tsv")

subprocess.run(
    ["hadoop", "fs", "-get", f"/tmp/{NAME}_QC2.tsv", f"../tmp/{NAME}_QC2.tsv"],
    check=True,
    shell=False,
)

# %%
