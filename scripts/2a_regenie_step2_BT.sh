#!/bin/sh

TRAIT="BT"
data_file_dir="/mnt/project/Data"

prompt="Enter GENE and PHENOTYPE for Step 2 (GENE burden files, PHENOTYPE .tsv file, and PHENOTYPE.${TRAIT}.LOCO/FIRTH files must exist):   "
read -p "$prompt" GENE PHENOTYPE 

run_regenie_step2="
sed -i "s|/home/dnanexus/out/out/|${data_file_dir}/step1/${PHENOTYPE}.${TRAIT}.LOCO/| ${PHENOTYPE}.${TRAIT}.step1_pred.list &&
sed -i "s|/home/dnanexus/out/out/|${data_file_dir}/step1/${PHENOTYPE}.${TRAIT}.LOCO/|' ${PHENOTYPE}.${TRAIT}.step1_firth.list && 
regenie \
  --step 2 \
  --bt \
  --bgen ${data_file_dir}/burden/${GENE}.bgen \
  --sample ${data_file_dir}/burden/${GENE}.sample \
  --ref-first \
  --phenoFile ${data_file_dir}/phenotypes/${PHENOTYPE}.${TRAIT}.final.tsv \
  --covarFile ${data_file_dir}/phenotypes/covariates.tsv \
  --firth --approx \
  --use-null-firth ${PHENOTYPE}.${TRAIT}.step1_firth.list \
  --pred ${PHENOTYPE}.${TRAIT}.step1_pred.list \
  --bsize 200 \
  --anno-file ${data_file_dir}/burden/${GENE}.annotations \
  --set-list ${data_file_dir}/burden/${GENE}.setlist \
  --mask-def ${data_file_dir}/burden/${GENE}.masks \
  --aaf-bins 0.01,0.05 \
  --check-burden-files \
  --write-mask-snplist \
  --out ${PHENOTYPE}.${TRAIT}.${GENE}.step2 \
  #--maxstep-null 1 \
  #--maxiter-null 25000 \
  --verbose
"

dx run swiss-army-knife -iin="/Data/step1/${PHENOTYPE}.${TRAIT}.LOCO/${PHENOTYPE}.${TRAIT}.step1_pred.list" \
   -iin="/Data/step1/${PHENOTYPE}.${TRAIT}.LOCO/${PHENOTYPE}.${TRAIT}.step1_firth.list" \
   -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/Data/results/${PHENOTYPE}.${TRAIT}.${GENE}" --yes;