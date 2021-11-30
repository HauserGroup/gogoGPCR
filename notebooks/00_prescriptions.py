# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.2
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
import pyspark

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
scripts = dataset["gp_scripts"]

field_names = [field.name for field in scripts.fields]
df = scripts.retrieve_fields(names=field_names, engine=dxdata.connect())

# %%
df.show(10, truncate=False)

# %%
df = spark.sql(
    "SELECT * FROM gp_scripts WHERE bnf_code IS NOT NULL AND LENGTH(bnf_code) = 14"
)
print(f"Rows with valid BNF codes: {df.count()}")
df.show(10, truncate=False)
