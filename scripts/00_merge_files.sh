#!/bin/sh

# Requirements: 
# Please refer to readme.md for more information about the requirements

# How to Run:
# Run this script using: sh 01-merge-files-dxfuse.sh on the command line

# Inputs
# - empty dummy.file

# Outputs
# - ukb_allChrs.bed - used as input for 1a/1b
# - ukb_allChrs.bim - used as input for 1a/1b
# - ukb_allChrs.fam - used as input for 1a/1b

#cmd to run (use as input with `-icmd=$run_merge`)
run_merge="cp /mnt/project/Bulk/Genotype\ Results/Genotype\ calls/ukb22418_c[1-9]* . ;\
        ls *.bed | sed -e 's/.bed//g'> files_to_merge.txt; \
        plink --merge-list files_to_merge.txt --make-bed\
        --autosome-xy --out ukb_allChrs;\
        rm files_to_merge.txt;\
        rm ukb22418_c*_b*_v2.*;"

dx run swiss-army-knife -iin="/Data/dummy.file" \
   -icmd="${run_merge}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/Data/step1" --brief --yes 