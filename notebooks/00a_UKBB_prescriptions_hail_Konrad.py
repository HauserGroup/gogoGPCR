# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.1
#   kernelspec:
#     display_name: Python [default]
#     language: python
#     name: python3
# ---

# %% [markdown]
# # UKBB Prescription Data

# %% [markdown]
# ## Load Packages

# %%
# Import packages
import os
import pandas as pd
import numpy as np
import hail as hl
import re

import datetime
from hail.plot import show
from pprint import pprint
import scipy
import statistics as st

hl.init()
hl.plot.output_notebook()

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
import pandas as pd
import re
import pickle

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

# %% [markdown]
# ## Import Data into Hail

# %%
"""
Field names are originally from Courtney's pandas code
"""

raw_data_loc = "gs://ukb31063/ukb31063.gp_scripts.20191008.txt"
# fs = [f"f{i}" for i in range(0, 8)]
# fields = ['eid','data_provider','issue_date','read_v2','bnf','dmd','drug_name','drug_quantity']
hl_presc = hl.import_table(raw_data_loc, no_header=False, delimiter="\t", impute=True)

drug_counts = hl_presc.group_by(hl_presc.drug_name).aggregate(counts=hl.agg.count())
drug_counts.show()

# %% [markdown]
# ## Export Drugs with >10k, >1k, and 100 Counts

# %%
# For the initial curation task, I restricted to the ~700 scripts
# that had >10,000 instances
top_drugs_10k = drug_counts.filter(drug_counts.counts > 10000)
top_drugs_10k.export("gs://gsarma/top_drugs_10k.tsv")

# Drugs that have >1000 scripts likely cover what we
# might want for future studies. An efficient way to bootstrap
# from the manually curated list of ~700 to the ~3500 in this
# larger list is a significant step.
top_drugs_1k = drug_counts.filter(drug_counts.counts > 1000)
top_drugs_1k.export("gs://gsarma/top_drugs_1k.tsv")

# %% [markdown]
# # Pre-curation Steps

# %%
"""
This is the first step of the curation process, simply split on the first non-word character
(most often a whitespace), and extract the first token. This is typically the generic name, 
but each drug needs to be looked at and potentially corrected.  

I found that this part was faster to do in pandas than with Hail, largely because of the
availability of questions on StackExchange. 
"""


def fill_possible_generic(dfname):
    """
    Use first token in raw prescription string as candidate generic name
    """
    dfname = dfname.drop(
        0
    )  # The first row is the empty string, i.e. there are a large number of empty prescriptions
    expr = "(^[\w]+)"
    newdf = pd.DataFrame()
    newdf["original_prescription"] = dfname["drug_name"]
    newdf["possible_generic"] = newdf["original_prescription"].str.extract(
        expr, expand=True
    )
    newdf["possible_generic"] = newdf["possible_generic"].str.lower()

    return newdf


def fill_delivery_system(dfname):
    """
    Extract delivery system for drug such as tablet, capsule, inhaler, etc.
    """
    med_type = []
    for values in dfname["original_prescription"]:
        # PO
        if re.search("tab|tabs|tablet", values, re.IGNORECASE):
            med_type.append("Tablet")
        elif re.search("cap|caps|capsule", values, re.IGNORECASE):
            med_type.append("Capsule")
        elif re.search("oral liquid", values, re.IGNORECASE):
            med_type.append("Oral Liquid")
        elif re.search("oral suspension", values, re.IGNORECASE):
            med_type.append("Oral Suspension")
        elif re.search("oral powder", values, re.IGNORECASE):
            med_type.append("Oral Powder")
        elif re.search("solution", values, re.IGNORECASE):
            med_type.append("Solution")
        elif re.search("sublingual spray", values, re.IGNORECASE):
            med_type.append("Sublingual Spray")
        elif re.search(
            "oral liquid|liquid peppermint|oral lyophilisates", values, re.IGNORECASE
        ):
            med_type.append("Oral Liquid")
        elif re.search("ispaghula husk", values, re.IGNORECASE):
            med_type.append("Food Supplement")

        # Nasal
        elif re.search("nasal cream", values, re.IGNORECASE):
            med_type.append("Nasal Cream")
        elif re.search("nasal spray", values, re.IGNORECASE):
            med_type.append("Nasal Spray")

        # Ear
        elif re.search("ear spray", values, re.IGNORECASE):
            med_type.append("Ear Spray")
        elif re.search("ear drops", values, re.IGNORECASE):
            med_type.append("Ear Drops")

        # Eye
        elif re.search("eye ointment", values, re.IGNORECASE):
            med_type.append("Eye Ointment")
        elif re.search("eye drops", values, re.IGNORECASE):
            med_type.append("Eye Drops")

        # Topical
        elif re.search("gel", values, re.IGNORECASE):
            med_type.append("Gel")
        elif re.search("ointment", values, re.IGNORECASE):
            med_type.append("Ointment")
        elif re.search("cream|emollient", values, re.IGNORECASE):
            med_type.append("Cream")
        elif re.search("shampoo", values, re.IGNORECASE):
            med_type.append("Shampoo")
        elif re.search("lotion|shower emollient", values, re.IGNORECASE):
            med_type.append("Lotion")
        elif re.search("scalp application", values, re.IGNORECASE):
            med_type.append("Cream")
        elif re.search("medicated nail lacquer", values, re.IGNORECASE):
            med_type.append("Medicated Nail Lacquer")
        elif re.search("patch", values, re.IGNORECASE):
            med_type.append("Patch")
        elif re.search("vaginal moisturiser", values, re.IGNORECASE):
            med_type.append("Vaginal Moisturizer")

        # Diagnostic
        elif re.search("testing strip", values, re.IGNORECASE):
            med_type.append("Testing Strip")

        # Other
        elif re.search(
            "inhaler|Evohaler|Accuhaler|Turbohaler|inhalation", values, re.IGNORECASE
        ):
            med_type.append("Inhaler")
        elif re.search("disc", values, re.IGNORECASE):
            med_type.append("Disc Inhaler")
        elif re.search("volumatic|AeroChamber", values, re.IGNORECASE):
            med_type.append("Inhaler Assist Device")
        elif re.search("pessary|pessaries", values, re.IGNORECASE):
            med_type.append("Pessary")
        elif re.search("needle|injection", values, re.IGNORECASE):
            med_type.append("Injection")
        elif re.search("vaccine", values, re.IGNORECASE):
            med_type.append("Vaccine")
        elif re.search("lancet", values, re.IGNORECASE):
            med_type.append("Lancet")

        else:
            med_type.append("NA")

    dfname["Delivery_System"] = med_type

    return dfname


def fill_dosages(dfname):
    """
    Extract dosages. Note that fill_delivery_system must be run BEFORE
    fill_dosages can be run, as the dosage patterns are
    dependent on the type of delivery system.
    """
    dosages = []
    single_dose_pat = "[0-9]+(\.)*[0-9]*(([a-zA-z]*)|(\%)*)"
    cream_pat = "[0-9]+(\.)*[0-9]*%"

    for values in dfname.iterrows():
        dtype = values[1]["Delivery_System"]
        presc = values[1]["original_prescription"]
        if (dtype == "Tablet") | (dtype == "Capsule") | (dtype == "Eye Drops"):
            if re.search(single_dose_pat, presc):
                dosages.append(re.search(single_dose_pat, presc).group(0))
            else:
                dosages.append("NA")
        elif dtype == "Cream":
            if re.search(cream_pat, presc):
                dosages.append(re.search(cream_pat, presc).group(0))
            else:
                dosages.append("NA")
        else:
            dosages.append("NA")

    dfname["Dosage"] = dosages

    return dfname


def pre_curation(dfname):
    """
    Prepare the UKBB GP data for curation by filling in possible
    generic names, delivery systems, and dosages.
    """
    return fill_dosages(fill_delivery_system(fill_possible_generic(dfname)))


# %%
td10k = pd.read_csv("/Users/gsarma/Dropbox/Broad/UKBB/top_drugs_10k.tsv", sep="\t")
td1k = pd.read_csv("/Users/gsarma/Dropbox/Broad/UKBB/top_drugs_1k.tsv", sep="\t")
curated = pd.read_csv(
    "/Users/gsarma/git/UKBB_prescriptions/UKBB_prescriptions_v2.csv", sep=","
)

td1ka = pre_curation(td1k)

# %%
pre_curation(td10k)

# %% [markdown]
# # Bootstrapping From Curated Entries

# %%
top_drugs = pd.read_csv("/Users/gsarma/Dropbox/Broad/UKBB/top_drugs_1k.tsv", sep="\t")
top_drugs = top_drugs.drop(0)
top_drugs = top_drugs.sort_values(by="counts", ascending=False)
last = len(top_drugs.index) - 1

# %%
all_dist = []
for n in range(round(last / 2)):
    base = top_drugs.iloc[
        last - n
    ].drug_name  # Start at the end and work to the half-way point
    temp_dist = []
    for k in range(round(last / 2)):
        l_dist = fuzz.token_sort_ratio(base, top_drugs.iloc[k].drug_name)
        # Start at the *beginning* and work to the half-way point
        temp_dist.append(l_dist)
    m = max(temp_dist)
    # First find all the indices where a maximum happens and take the minimum
    all_dist.append(min([i for i, j in enumerate(temp_dist) if j == m]))

# %%
plt.hist(all_dist, bins=100)

# %%
check_quality = []
for n in range(0, 1700, 100):
    check_quality.append(
        [top_drugs.iloc[last - n].drug_name, top_drugs.iloc[all_dist[n]].drug_name]
    )

# %%
qc_df = pd.DataFrame(check_quality)
qc_df


# %%
def find_best_match(raw_script, reference):
    """
    Find the best match for a prescription string starting
    from a curated reference list.
    TODO: option to use either generic name or full prescription string
    """
    all_dist = []
    num_ref = len(reference) - 1
    token = raw_script.split(" ", 1)
    for k in range(num_ref):
        l_dist = fuzz.token_sort_ratio(raw_script, reference.iloc[k]["Generic_Name"])
        all_dist.append(l_dist)
    m = max(all_dist)
    # First find all the indices where a maximum happens and take the minimum
    index = min([i for i, j in enumerate(all_dist) if j == m])
    return [
        reference.iloc[index]["Generic_Name"],
        reference.iloc[index]["Drug_Category_and_Indication"],
    ]


def bootstrap_curation(new_df, reference):
    """
    Take in an uncurated list of prescriptions and bootstrap
    """
    candidate_generic = []
    candidate_category = []
    for values in new_df["possible_generic"]:
        best_match = find_best_match(values, curated)
        candidate_generic.append(best_match[0])
        candidate_category.append(best_match[1])

    return [candidate_generic, candidate_category]


# %%
print(bootstrap_curation.__doc__)

# %%
short_df = td1ka.sample(n=20)
new_columns = bootstrap_curation(short_df, curated)
short_df["possible_generic"] = new_columns[0]
short_df["Drug_Category_and_Indication"] = new_columns[1]

# %%
short_df

# %% [markdown]
# # Some Basic Data Exploration

# %%
plt.hist(top_drugs.drug_name.str.len())
plt.xlabel("Prescription String Legnth")
plt.ylabel("Frequency")
plt.legend()
plt.show()

# %%
td1ka

# %%
