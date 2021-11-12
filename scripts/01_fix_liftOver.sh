#!/bin/sh

# The GCRh37 to GRCh38 chain file contains alternate contigs which are imcompatible with Hail, this script fixes it in an ugly way

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

python fix_liftOver.py

rm hg19ToHg38.over.chain.gz