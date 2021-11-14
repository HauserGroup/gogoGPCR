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
import subprocess
from distutils.version import LooseVersion

import dxdata
import dxpy
import pyspark
from pyspark.sql import functions as F
from pyspark.sql.types import IntegerType

# %%
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)

dispensed_database_name = dxpy.find_one_data_object(
    classname="database", name="app*", folder="/", name_mode="glob", describe=True
)["describe"]["name"]
dispensed_dataset_id = dxpy.find_one_data_object(
    typename="Dataset", name="app*.dataset", folder="/", name_mode="glob"
)["id"]

dataset = dxdata.load_dataset(id=dispensed_dataset_id)
participant = dataset["participant"]


def fields_for_id(field_id):

    field_id = str(field_id)
    fields = participant.find_fields(
        name_regex=r"^p{}(_i\d+)?(_a\d+)?$".format(field_id)
    )

    return sorted(fields, key=lambda f: LooseVersion(f.name))


# %%
fields = [
    "21022",
    "22001",
    "22009",
]
field_names = [fields_for_id(id) for id in fields]
field_names = ["eid"] + [field.name for fields in field_names for field in fields]

pcs = {f"p22009_a{i}": f"PC{i}" for i in range(1, 21)}
covs = ["FID", "IID", "AGE", "AGE2", "AGESEX", "AGE2SEX"] + list(pcs.values())

# %%
df = participant.retrieve_fields(
    names=field_names, engine=dxdata.connect(), coding_values="raw"
)

df = df.na.drop(how="any")

df = (
    df.select([F.col(c).alias(pcs.get(c, c)) for c in df.columns])
    .withColumn("FID", F.col("eid"))
    .withColumn("IID", F.col("eid"))
    .withColumn("AGE", F.col("p21022").cast(IntegerType()))
    .withColumn("AGE2", (F.col("p21022") ** 2).cast(IntegerType()))
    .withColumn("AGESEX", (F.col("p21022") * F.col("p22001")).cast(IntegerType()))
    .withColumn(
        "AGE2SEX", ((F.col("p21022") ** 2) * F.col("p22001")).cast(IntegerType())
    )
    .select(*covs)
)

df.show(5, truncate=False)

# %%
df.coalesce(1).write.csv(
    "/tmp/covariates.tsv",
    sep="\t",
    header=True,
)

# %%
subprocess.run(
    ["hadoop", "fs", "-getmerge", "/tmp/covariates.tsv", "../tmp/covariates.tsv"],
    check=True,
    shell=False,
)
subprocess.run(
    ["dx", "upload", "../tmp/covariates.tsv", "--path", "/Data/phenotypes/"],
    check=True,
    shell=False,
)
