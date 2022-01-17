# Notebooks

Notebooks should generally be run in order. Notebooks starting with `0*_` need only be run once per configuration and can be re-used for analyses.

0. `00_prescriptions.py`: Fetch all prescription records for EIDs in `Participant_table.csv`. **Not run as part of pipeline**
1. `01_sample_hard_filter.py`: Generate list of samples to filter out based quality measures and  relatedness. Requires `ukb_rel.dat` file found buried somewhere within UKB Genotyping Quality folders. Threshold can be set in `config.toml`. Can be re-used across analyses.
2. `02_covariates.py`: Generate standard covariates as per *Mbatchou et al., 2021* for all 500k samples. Columns: FID, IID, sex, age, age*sex, age^2, age^2*sex, and 20 first PCs. Can be re-used across analyses.
3. `03a_psychiatric_phenotypes_BT.py` Generate PHESANT-friendly phenotype TSV file for psychiatric ICD10 diagnoses with more than 500 cases. Phenotype scripts are bespoke. 
4. `03b_metabolic_phenotypes_QT.py` Quantitative trait (QT) metabolic phenotypes for 200k WES samples.
5. `03c_metabolic_phenotypes_BT.py` Binary trait (BT) metabolic phenotypes.
6. `04_fix_phesant_output.py` fixes the output from `scripts/04_run_PHESANT.sh`.
7. `1_qc1_matrixtables.py` performs initial QC on raw VCF files. See notebook for specific steps.
8. `2a/b_annotate_*.py` annotates variants in MatrixTable from QC1 for burden testing. This step is necessary for running regenie and is bespoke to your particular analysis. 
9. `3_qc2_matrixtables.py` Performs further QC on variants, samples, and genotyping quality. Generates .annotations, .setlist, .bgen, and .sample file necessary for regenie. See notebook for detailed explanation of steps.
10. `4_forest_plots.py` visualise analysis results output from regenie step 2, e.g. `2a_regenie_step2_BT.sh` as Forest plots. See cited article from example of plots.
