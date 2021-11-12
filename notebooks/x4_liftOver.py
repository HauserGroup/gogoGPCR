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
import pyspark
import dxpy
import hail as hl
import gzip

# %%
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)
hl.init(sc=sc, default_reference="GRCh38")

# %%
rg37 = hl.get_reference("GRCh37")
rg38 = hl.get_reference("GRCh38")
rg37.add_liftover(
    "file:" + "/opt/notebooks/gogoGPCR/tmp/hg19ToHg38.over.filtered.chain.gz", rg38
)

# %%
mt = hl.import_vcf(
    "file:"
    + "/mnt/project/Bulk/Exome\ sequences/Population\ level\ exome\ OQFE\ variants\,\ pVCF\ format/ukb23156_c13_b7_v1.vcf.gz",
    force_bgz=True,
    array_elements_required=False,
    drop_samples=True,
)

# %%
mt.checkpoint("/tmp", overwrite=True)

# %%
mt.show()


# %%
def rreplace(s, old, new, occurrence):
    li = s.rsplit(old, occurrence)
    return new.join(li)


writing = True

with gzip.open(
    "/opt/notebooks/gogoGPCR/tmp/hg19ToHg38.over.filtered.chain.gz", "wt"
) as new:
    for f in gzip.open("/opt/notebooks/gogoGPCR/tmp/hg19ToHg38.over.chain.gz", "rt"):
        if f.startswith("chain") and (("_" in f) or ("M" in f)):
            writing = False
        if f.startswith("chain") and (("_" not in f) and ("M" not in f)):
            writing = True
            f = f.replace("chr", "", 1)

        if writing:
            new.write(f)


# %%
# str.replace?

# %%
mt.count()

# %%
