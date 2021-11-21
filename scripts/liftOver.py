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
import pyspark
import dxpy
import hail as hl
import subprocess

# %%
sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)
hl.init(sc=sc, log="../tmp/hail_logs/liftover_log.txt")

# %%
# subprocess.run(["hadoop,", "fs", "-rm", "-r", "/tmp",], check = True, shell = False)
# subprocess.run(["hadoop", "fs", "-mkdir", "/tmp"], check = True, shell = False)
# subprocess.run(["hadoop", "fs", "-put", "../data/misc/hg19ToHg38.over.filtered.chain.gz", "/tmp/" ], check = True, shell = False)

# %%
rg37 = hl.get_reference("GRCh37")
rg38 = hl.get_reference("GRCh38")
rg37.add_liftover("/tmp/hg19ToHg38.over.filtered.chain.gz", rg38)

# %%
allChrs = "file:" + "/mnt/project/Data/step1/ukb_allChrs"

mt = hl.import_plink(
    bed=allChrs + ".bed",
    bim=allChrs + ".bim",
    fam=allChrs + ".fam",
    reference_genome="GRCh37",
    min_partitions=800,
    missing="-9",
)


# %%
mt.checkpoint("/tmp/cp.mt", overwrite=True)

# %%
mt = mt.annotate_rows(
    new_locus=hl.liftover(mt.locus, rg38, include_strand=True),
    old_locus=mt.locus,
)
mt = mt.filter_rows(
    hl.is_defined(mt.new_locus) & ~mt.new_locus.is_negative_strand
)
print(f"New number of loci: {mt.count()}")
mt = mt.key_rows_by(locus=mt.new_locus.result, alleles=mt.alleles)

# %%
hl.export_plink(
    mt, "/tmp/ukb_allChrs.GRCh38", fam_id=mt.s, is_female=mt.is_female, pheno=-9
)

# %%
# subprocess.run(["hadoop", "fs", "-get", "/tmp/ukb_allChrs.GRCh38.*"], check = True, shell = False)
# subprocess.run(["dx", "upload", "ukb_allChrs.GRCh38.*", "--path", "/Data/step1/"], check = True, shell = False)
