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

# %%
import re
import subprocess

import dxdata
import dxpy
import pandas as pd
import pyspark
from pyspark.sql import functions as F

# %%
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

dispensed_database_name = dxpy.find_one_data_object(
    classname="database", name="app*", folder="/", name_mode="glob", describe=True
)["describe"]["name"]
dispensed_dataset_id = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob"
)["id"]

spark.sql("USE " + dispensed_database_name)

dataset = dxdata.load_dataset(id=dispensed_dataset_id)
participant = dataset["participant"]


def get_columns_to_keep(df, threshold=200) -> list:
    """
    This function drops all columns which contain more null values than threshold
    :param df: A PySpark DataFrame
    """
    null_counts = (
        df.select([F.count(F.when(~F.col(c).isNull(), c)).alias(c) for c in df.columns])
        .collect()[0]
        .asDict()
    )
    to_keep = [k for k, v in null_counts.items() if v > threshold]
    # df = df.select(*to_keep)

    return to_keep


def new_names(s: str) -> str:
    """
    Fixes a column header for PHESANT use
    """
    s = s.replace("p", "x").replace("i", "")

    if bool(re.search("_\d$", s)):
        s += "_0"
    else:
        s += "_0_0"
    return s


# %%
first_occurences = list(
    participant.find_fields(
        lambda f: bool(re.match("^Date [F]\d{2} first reported", f.title))
    )
)

age_sex = "|".join(["31", "21022"])
age_sex_fields = list(
    participant.find_fields(lambda f: bool(re.match(f"^p({age_sex})$", f.name)))
)

psych = "|".join(
    [
        "2090",
        "2100",
        "20126",
    ]
)

psych_fields = list(
    participant.find_fields(lambda f: bool(re.match(f"^p({psych})\D", f.name)))
)

field_names = (
    ["eid"]
    + [f.name for f in age_sex_fields]
    + [f.name for f in first_occurences]
    + [f.name for f in psych_fields]
)
df = participant.retrieve_fields(names=field_names, engine=dxdata.connect())

# %%
min_cases = 500

print(f"Number of columns: {len(df.columns)}")
to_keep = get_columns_to_keep(df, min_cases)

to_keep.insert(1, to_keep[11])
to_keep.insert(2, to_keep[11])
to_keep.pop(12)
to_keep.pop(12)


df = df.select(*to_keep)
print(f"Number of columns with at least {min_cases} cases: {len(df.columns)}")

# %%
new_names = ["xeid"] + [new_names(s) for s in df.columns[1:]]
print(new_names[:10])

# %%
df = df.toDF(*new_names)

# %%
df.write.csv("/tmp/phenos.csv", sep=",", header=True)

# %%
subprocess.run(
    ["hadoop", "fs", "-rm", "/tmp/phenos.csv/_SUCCESS"], check=True, shell=False
)
subprocess.run(
    ["hadoop", "fs", "-get", "/tmp/phenos.csv", "phenos.csv"], check=True, shell=False
)

# %%
# !sed -e '3,${/^xeid/d' -e '}' phenos.csv/part* > phenos.filtered.csv
