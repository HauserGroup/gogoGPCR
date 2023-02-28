#!/bin/sh

TRAIT="BT"
data_file_dir="/mnt/project/Data"

prompt="Enter GENE and PHENOTYPE for Step 2 (GENE burden files, PHENOTYPE .tsv file, and PHENOTYPE.${TRAIT}.LOCO/FIRTH files must exist):   "
read -p "$prompt" GENE PHENOTYPE 

cp "${data_file_dir}/step1/${PHENOTYPE}.${TRAIT}.LOCO/${PHENOTYPE}.${TRAIT}.step1_firth.list" .
cp "${data_file_dir}/step1/${PHENOTYPE}.${TRAIT}.LOCO/${PHENOTYPE}.${TRAIT}.step1_pred.list" .

sed -i "s|/home/dnanexus/out/out/|${data_file_dir}/step1/${PHENOTYPE}.${TRAIT}.LOCO/|" "${PHENOTYPE}.${TRAIT}.step1_pred.list" 
sed -i "s|/home/dnanexus/out/out/|${data_file_dir}/step1/${PHENOTYPE}.${TRAIT}.LOCO/|" "${PHENOTYPE}.${TRAIT}.step1_firth.list"

mkdir -p "${PHENOTYPE}.${TRAIT}.${GENE}"

regenie \
  --step 2 \
  --bt \
  --bgen "${data_file_dir}/burden/${GENE}.bgen" \
  --sample "${data_file_dir}/burden/${GENE}.sample" \
  --ref-first \
  --phenoFile "${data_file_dir}/phenotypes/${PHENOTYPE}.${TRAIT}.final.tsv" \
  --covarFile "${data_file_dir}/phenotypes/covariates.tsv" \
  --firth \
  --use-null-firth "${PHENOTYPE}.${TRAIT}.step1_firth.list" \
  --pred "${PHENOTYPE}.${TRAIT}.step1_pred.list" \
  --bsize 200 \
  --anno-file "${data_file_dir}/burden/${GENE}.annotations" \
  --set-list "${data_file_dir}/burden/${GENE}.setlist" \
  --mask-def "${data_file_dir}/burden/${GENE}.masks" \
  --aaf-bins 0.0195 \
  --check-burden-files \
  --write-mask-snplist \
  --maxstep-null 2 \
  --maxiter-null 20000 \
  --out "${PHENOTYPE}.${TRAIT}.${GENE}/${PHENOTYPE}.${TRAIT}.${GENE}.step2" \
  --vc-tests acato-full \
  --vc-MACthr 10 \
  --verbose