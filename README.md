# gogoGPCR

gogoGPCR is a framework for performing burden testing on UK Biobank Whole-Exome Sequencing (WES) data. This repo contains a series of notebooks for pre-processing and quality-controlling UKB VCF files, mainly with [Hail](https://hail.is), a Dockerfile for pre-processing phenotypes with [PHESANT](https://github.com/MRCIEU/PHESANT), and a series of Shell scripts for performing burden testing, with [regenie](https://rgcgithub.github.io/regenie/).

For information on running individual notebooks, see `notebooks/README.md`. Likewise, see `scripts/README.md` for running individual scripts.

**WARNING: gogoGPCR will probably (99% confidence) not work out of the box and requires tweaking for your particular setup. In the very least, it may serve as inspiration for your own analyses.**

## Usage
Notebooks and scripts should generally be run in numerical order. Notebook/scripts starting with `0*_` e.g. `notebooks/02_covariates.py`, need only be run once and the generated file can be re-used for further analyses. Files wigth `*a,b,c_` e.g. `03c_metabolic_phenotypes_BT.py` represent different versions of the same step. 

## Structure

```
├── README.md
├── config.toml                     # Majority of configuration
├── notebooks                       # Pre-processing and QC Jupyter notebooks
│   ├── 00_prescriptions.py         
│   ├── 01_sample_hard_filter.py    
│   ├── ...
│   ├── 4_forest_plots.py
│   └── README.md
├── docker                          # Dockerfile for running PHESANT
│   └── phesant  
│       ├── install_packages.R  
│       └── Dockerfile
├── scripts                         # Shell scripts for running regenie step 1 and 2
│   ├── 00_merge_files.sh
│   ├── ...
│   ├── 1a_regenie_step1_BT.sh
│   ├── 2a_regenie_step2_BT.sh
│   └── README.md    
├── src                             # QC and utility functions
│   ├── matrixtables.py
│   ├── ...
│   └── utils.py  
└── data
    ├── examples                    # Example regenie output
    │   └── *.regenie
    └── misc                        # Miscellaneous data files for mapping purposes
        ├── Data_Dictionary_Custom.csv
        ├── ...
        └── xgen_plus_spikein.b38.bed

```

## Prefer .ipynb notebooks?
Run `jupytext --from py:percent --to notebook notebooks/*.py`


## Citation
<div class="csl-entry">van der Velden, W. J. C., Lindquist, P., Madsen, J. S., Stassen, R. H. M. J., Wewer Albrechtsen, N. J., Holst, J. J., Hauser, A. S., &#38; Rosenkilde, M. M. (2021). Molecular and in vivo phenotyping of missense variants of the human glucagon receptor. <i>Journal of Biological Chemistry</i>, <i>0</i>(0), 101413. https://doi.org/10.1016/j.jbc.2021.101413</div>