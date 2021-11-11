#!/bin/sh

# This script takes roughly 1 hour on 16x8 (cores x nodes) and converts the merged PLINK files to GRCh38.
# Should /Data/step1/ukb_allChrs.GRCh38.(bed, bim, fam) already exist. Don't re-run.
# Don't be scared if Spark stages fail once or twice. Be scared if they fail 4 times. 

hadoop fs -rm -r /tmp
hadoop fs -mkdir /tmp
hadoop fs -put ../data/misc/hg19ToHg38.over.filtered.chain.gz /tmp/

python liftOver.py

echo "copying files from HDFS to local"

hadoop fs -get /tmp/ukb_allChrs.GRCh38.*

echo "uploading files from local to project storage"

dx upload ukb_allChrs.GRCh38.* --path /Data/step1/