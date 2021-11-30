#!/bin/sh

gene="DRD2"
phenotype="psychiatric"
data_file_dir="/mnt/project/Data"

run_regenie_step2="
sed -i 's|/home/dnanexus/out/out/|${data_file_dir}/step1/${phenotype}.BT/|' ${phenotype}_step1_BT_pred.list &&
sed -i 's|/home/dnanexus/out/out/|${data_file_dir}/step1/${phenotype}.BT/|' ${phenotype}_step1_BT_firth.list && 
regenie \
  --step 2 \
  --bt \
  --bgen ${data_file_dir}/burden/${gene}.bgen \
  --sample ${data_file_dir}/burden/${gene}.sample \
  --ref-first \
  --phenoFile ${data_file_dir}/phenotypes/${phenotype}.tsv \
  --covarFile ${data_file_dir}/phenotypes/covariates.tsv \
  --firth --approx \
  --use-null-firth ${phenotype}_step1_BT_firth.list \
  --pred ${phenotype}_step1_BT_pred.list \
  --bsize 200 \
  --anno-file ${data_file_dir}/burden/${gene}.annotations \
  --set-list ${data_file_dir}/burden/${gene}.setlist \
  --mask-def ${data_file_dir}/burden/${gene}.masks \
  --check-burden-files \
  --out step2_BT_${phenotype}_${gene} \
  #--maxstep-null 1 \
  #--maxiter-null 25000 \
  --write-samples \
  --print-pheno \
  --af-cc \
  --verbose
"

dx run swiss-army-knife -iin="/Data/step1/${phenotype}.BT/psychiatric_step1_BT_pred.list" \
   -iin="/Data/step1/${phenotype}.BT/psychiatric_step1_BT_firth.list" \
   -icmd="${run_regenie_step2}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/Data/results" --yes;