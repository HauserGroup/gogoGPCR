import pyspark
import dxpy
import hail as hl
import gzip
import subprocess

sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)
hl.init(sc=sc, default_reference='GRCh38')

#do wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz subprocess.run(["dx", "upload", bgen_file, sample_file, ANNOTATIONS_FILE, SETLIST_FILE, "--path", "/data/burden/"], check = True, shell = False)

#rg37 = hl.get_reference('GRCh37')  
#rg38 = hl.get_reference('GRCh38')  
#rg37.add_liftover("file:" + "/opt/notebooks/gogoGPCR/tmp/hg19ToHg38.over.filtered.chain.gz", rg38)

writing = True

with gzip.open("/opt/notebooks/gogoGPCR/tmp/hg19ToHg38.over.filtered.chain.gz", "wt") as new:
    for f in gzip.open("/opt/notebooks/gogoGPCR/tmp/hg19ToHg38.over.chain.gz", "rt"):
        if f.startswith("chain") and (("_" in f) or ("M" in f)):
            writing = False
        if f.startswith("chain") and (("_" not in f) and ("M" not in f)):
            writing = True
            f = f.replace("chr", "", 1)
            
        if writing:
            new.write(f)
# do dx upload