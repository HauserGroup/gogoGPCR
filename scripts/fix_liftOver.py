import gzip
import subprocess

# subprocess.run(["wget", "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"], check = True, shell = False)

writing = True

with gzip.open("../data/misc/hg19ToHg38.over.filtered.chain.gz", "wt") as new:
    for f in gzip.open("hg19ToHg38.over.chain.gz", "rt"):
        if f.startswith("chain") and (("_" in f) or ("M" in f)):
            writing = False
        if f.startswith("chain") and (("_" not in f) and ("M" not in f)):
            writing = True
            f = f.replace("chr", "", 1)

        if writing:
            new.write(f)

# subprocess.run(["rm", "hg19ToHg38.over.chain.gz"], check = True, shell = False)