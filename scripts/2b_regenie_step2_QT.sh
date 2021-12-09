#!/bin/sh

TRAIT="QT"
data_file_dir="/mnt/project/Data"

prompt="Enter GENE and PHENOTYPE for Step 2 (GENE burden files, PHENOTYPE .tsv file, and PHENOTYPE.${TRAIT}.LOCO files must exist)"
read -p "$prompt" GENE PHENOTYPE 

run_regenie_step2="
#sed -i 's|/home/dnanexus/out/out/|${data_file_dir}/step1/${PHENOTYPE}.${TRAIT}.LOCO/|' ${PHENOTYPE}.${TRAIT}.step1_pred.list &&
regenie \
  --step 2 \
  --qt \
  --bgen ${data_file_dir}/burden/${GENE}.bgen \
  --sample ${data_file_dir}/burden/${GENE}.sample \
  --ref-first \
  --phenoFile ${data_file_dir}/phenotypes/${PHENOTYPE}.${TRAIT}.final.tsv \
  --covarFile ${data_file_dir}/phenotypes/covariates.tsv \
  --firth --approx \
  --pred ${PHENOTYPE}.${TRAIT}.step1_pred.list \
  --anno-file ${data_file_dir}/burden/${GENE}.annotations \
  --set-list ${data_file_dir}/burden/${GENE}.setlist \
  --mask-def ${data_file_dir}/burden/${GENE}.masks \
  --check-burden-files \
  --out ${PHENOTYPE}.${TRAIT}.${GENE}.step2 \
  #--maxstep-null 2 \
  #--maxiter-null 15000 \
  --write-samples \
  --print-pheno \
  --verbose
"

dx run swiss-army-knife -iin="/Data/step1/${PHENOTYPE}.${TRAIT}/${PHENOTYPE}.${TRAIT}.step1_pred.list" \
   -icmd="${run_regenie_step2}" --tag="Step2" --instance-type "mem1_ssd1_v2_x16"\
   --destination="/Data/results/${PHENOTYPE}.${TRAIT}.${GENE}" --yes;