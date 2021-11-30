#!/bin/sh

pheno_file_dir="/mnt/project/Data/phenotypes"
phenotype="psychiatric"

run_phesant="mkdir -p /tmp/${phenotype}.phesant && pwd && \
    cd /PHESANT/WAS && \
    sed -e '1s|^|setwd(\"/PHESANT/WAS\")\n|' phenomescan.R && \
    Rscript phenomeScan.r --phenofile='${pheno_file_dir}/${phenotype}.raw.csv' --variablelistfile='../variable-info/outcome-info.tsv' --datacodingfile='../variable-info/data-coding-ordinal-info.txt' --resDir='/tmp/${phenotype}.phesant' --userId='xeid'  --save"
        
dx run swiss-army-knife -iin="/Data/dummy.file" \
   -icmd="${run_phesant}" --tag="PHESANT" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/Data/phenotypes" \
   -iimage="jsture/phesant" --yes;