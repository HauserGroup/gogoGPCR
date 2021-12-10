#!/bin/sh

# Requirements: 
# 03_plink_qc.sh must have been run

# How to Run:
# Run this shell script using: 
#   sh partD-step1-qc-regenie.sh 
# on the command line on your own machine

# Inputs:
# Note that you can adjust the output directory by setting the data_file_dir variable
# - /Data/diabetes_wes_200k.phe - from part A (please refer to notebook & slides)
# - /Data/WES_array_snps_qc_pass.snplist - from part C
# - /Data/ukb22418_c1_22_merged.bed - from part B
# - /Data/ukb22418_c1_22_merged.bed - from part B
# - /Data/ukb22418_c1_22_merged.bed - from part B

# Outputs:
# - /Data/diabetes_results_1.loco.gz - Leave One Chromosome Out results (used in part F)
# - /Data/diabetes_results_pred.list - List of files generated this step (used in part F)
# - /Data/diabetes_results.log


#output directory - this should also be where the files in 02-step1-qc-filter.sh end up
TRAIT="BT"
step1_file_dir="/mnt/project/Data/step1"
pheno_file_dir="/mnt/project/Data/phenotypes"
data_file_dir="Data/step1"

prompt="Enter phenotype for Step 1 (PHENOTYPE.${TRAIT}.final.tsv must exist and contain only binary phenotypes):   "
read -p "$prompt" PHENOTYPE 

run_regenie_step1="
regenie \
 --step 1 --bt \
 --out ${PHENOTYPE}.${TRAIT}.step1 \
 --bed ukb_allChrs.GRCh38 \
 --phenoFile ${pheno_file_dir}/${PHENOTYPE}.${TRAIT}.final.tsv --covarFile ${pheno_file_dir}/covariates.tsv \
 --extract ${step1_file_dir}/WES_qc_pass.snplist --keep ${step1_file_dir}/WES_qc_pass.id \
 --bsize 1000 \
 --write-null-firth \
 --lowmem --lowmem-prefix tmp_preds \
 --verbose --threads 16
"

dx run swiss-army-knife -iin="${data_file_dir}/ukb_allChrs.GRCh38.bed" -iin="${data_file_dir}/ukb_allChrs.GRCh38.bim" -iin="${data_file_dir}/ukb_allChrs.GRCh38.fam" \
   -icmd="${run_regenie_step1}" \
   --tag="Step1" --instance-type "mem1_ssd1_v2_x16" --destination="/Data/step1/${PHENOTYPE}.${TRAIT}.LOCO" --brief --yes