#!/bin/sh

gene=$1
phenotype=$2

data_file_dir="/mnt/project/Data"

run_regenie_step2="
sed -i 's|/home/dnanexus/out/out/|${data_file_dir}/step1/${phenotype}.QT/|' ${phenotype}_step1_QT_pred.list &&
sed -i 's|/home/dnanexus/out/out/|${data_file_dir}/step1/${phenotype}.QT/|' ${phenotype}_step1_QT_firth.list && 
regenie \
  --step 2 \
  --qt \
  --bgen ${data_file_dir}/burden/${gene}.bgen \
  --sample ${data_file_dir}/burden/${gene}.sample \
  --ref-first \
  --phenoFile ${data_file_dir}/phenotypes/${phenotype}.tsv \
  --covarFile ${data_file_dir}/phenotypes/covariates.tsv \
  --phenoColList '48,95,102,3148,21001,21002,23099,23127,23225,23231,23236,30750,30760,30780,30870,diastolic,systolic,pulse' \
  --firth --approx \
  --pred ${phenotype}_step1_QT_pred.list \
  --anno-file ${data_file_dir}/burden/${gene}.annotations \
  --set-list ${data_file_dir}/burden/${gene}.setlist \
  --mask-def ${data_file_dir}/burden/${gene}.masks \
  --check-burden-files \
  --out step2_QT_${phenotype}_${gene} \
  #--maxstep-null 2 \
  #--maxiter-null 15000 \
  --write-samples \
  --print-pheno \
  --af-cc \
  --verbose
"

dx run swiss-army-knife -iin="/Data/step1/${phenotype}.QT/${phenotype}_step1_QT_pred.list" \
   -iin="/Data/step1/${phenotype}.QT/${phenotype}_step1_QT_firth.list" \
   -icmd="${run_regenie_step2}" --tag="Step1" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/Data/results/" --yes;