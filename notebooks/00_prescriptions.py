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

# %%
# IMPORTS
import re

import dxdata
import dxpy
import pyspark
from pprint import pprint

import tomli
import hail as hl
from pathlib import Path
from pyspark.sql.types import StringType
import subprocess

# SET LOGFILE

with open("../config.toml", "rb") as f:
    conf = tomli.load(f)

    LOG_FILE = Path(conf["IMPORT"]["LOG_DIR"], f"prescriptions.log").resolve().__str__()


# INIT SPARK
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

# INIT HAIL
hl.init(sc=sc, default_reference="GRCh38", log=LOG_FILE)

# DISPENSE DATASET
dispensed_database_name = dxpy.find_one_data_object(
    classname="database", name="app*", folder="/", name_mode="glob", describe=True
)["describe"]["name"]
dispensed_dataset_id = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob"
)["id"]

spark.sql("USE " + dispensed_database_name)

dataset = dxdata.load_dataset(id=dispensed_dataset_id)
scripts = dataset["gp_scripts"]

field_names = [field.name for field in scripts.fields]
df = scripts.retrieve_fields(names=field_names, engine=dxdata.connect())

# MAKE TMPDIR

Path("../tmp").resolve().mkdir(parents=True, exist_ok=True)

# %%
df = df.withColumn("issue_date", df.issue_date.cast(StringType()))

# %%
ht = hl.Table.from_spark(df, key=["eid"])

# %%
ids = hl.import_table(
    "file://" + "/mnt/project/Data/schizophrenia/Participant_table.csv",
    delimiter=",",
    impute=True,
    key="Participant ID",
    types={"Participant ID": "str"},
)

# %%
ht = ht.semi_join(ids)
ht = ht.annotate(**ids[ht.eid])
ht.show()

# %%
total = ht.count()
patients = ht.group_by("eid").aggregate(eids=hl.agg.count()).count()
pprint(f"Total number of prescriptions: {total} from {patients} unique patients")

# %%
write_path = "/tmp/schizophrenia_prescriptions.tsv.bgz"
ht.export(write_path)

# %%
subprocess.run(
    ["hadoop", "fs", "-get", write_path, f"..{write_path}"], check=True, shell=False
)

# %%
subprocess.run(
    [
        "dx",
        "upload",
        f"..{write_path}",
        "--path",
        "Data/schizophrenia/schizophrenia_prescriptions.tsv.bgz",
    ],
    check=True,
    shell=False,
)
