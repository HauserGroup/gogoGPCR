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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
import os
import sys
from pathlib import Path

import dxdata
import dxpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyspark
from matplotlib.lines import Line2D

module_path = Path("..").resolve().__str__()
if module_path not in sys.path:
    sys.path.append(module_path)

import src.results as sr


# %%
# Flags
GENE = "DRD2"
TRAIT = "QT"

# %%
# sc = pyspark.SparkContext()
# spark = pyspark.sql.SparkSession(sc)

# %%
# Results files
files = [
    f"file:/mnt/project/data/results/{GENE}/{file}"
    for file in os.listdir(f"/mnt/project/data/results/{GENE}")
]

# Field codings
ukb_coding = pd.read_csv(
    "/opt/notebooks/gogoGPCR/data/misc/Data_Dictionary_Showcase.csv",
    error_bad_lines=False,
    warn_bad_lines=False,
    quotechar='"',
    usecols=["FieldID", "Field"],
)

custom_coding = pd.read_csv(
    "/opt/notebooks/gogoGPCR/data/misc/Data_Dictionary_Custom.csv",
    error_bad_lines=False,
    warn_bad_lines=False,
    quotechar='"',
)

# %%
# Load raw DF
df_raw = pd.read_csv(files[0], delimiter=" ", header="infer", comment="#").assign(
    SOURCE=os.path.basename(files[0])
)
df_raw = pd.concat(
    [df_raw]
    + [
        pd.read_csv(fp, delimiter=" ", comment="#").assign(SOURCE=os.path.basename(fp))
        for fp in files[1:]
    ],
    axis=0,
)

# %%
# Fix common fields
df = df_raw
df.loc[:, "GENE"] = df.ID.apply(lambda x: x.split(".")[0])
df.loc[:, "MASK"] = df.ALLELE1.apply(lambda x: x.split(".", maxsplit=2)[0])
df.loc[:, "AAF"] = df.ALLELE1.apply(lambda x: x.split(".", maxsplit=1)[-1])
# df.loc[:,"MASK2"] = df.MASK + "." + df.AAF
df.loc[:, "FILE"] = df.SOURCE.apply(
    lambda x: x.split("/")[-1].split("\.")[0].split("_", maxsplit=4)[2:]
)
df.loc[:, "TRAIT"] = df.FILE.apply(lambda x: x[0])
df.loc[:, "PHENO"] = df.FILE.apply(lambda x: x[-1].split(".")[0])
df = df.drop(["ID", "ALLELE0", "ALLELE1", "EXTRA", "SOURCE", "FILE", "TEST"], axis=1)

# Sanity check
df.head()

# %%
# Filters
bt = df.TRAIT == "BT"
qt = df.TRAIT == "QT"

# Fix Binary Traits
df.loc[bt, "OR"] = np.exp(df.loc[bt, "BETA"])
df.loc[bt, "OR_up"] = np.exp(df.loc[bt, "BETA"] + df.loc[bt, "SE"])
df.loc[bt, "OR_low"] = np.exp(df.loc[bt, "BETA"] - df.loc[bt, "SE"])
df.loc[bt, "OR_up_lim"] = df.loc[bt, "OR_up"] - df.loc[bt, "OR"]
df.loc[bt, "OR_low_lim"] = df.loc[bt, "OR"] - df.loc[bt, "OR_low"]

# Fix Quantitative Traits
df.loc[qt, "BETA_up_lim"] = df.loc[qt, "BETA"] + df.loc[qt, "SE"]
df.loc[qt, "BETA_low_lim"] = df.loc[qt, "BETA"] - df.loc[qt, "SE"]

# Final fixes
df.loc[:, "Phenotype"] = df.PHENO.apply(
    lambda x: sr.pheno_search(x, ukb_coding, custom_coding).replace('"', "").strip()
)
df.loc[:, "pval"] = np.power(10, -df["LOG10P"])
df.loc[:, "pval_stars"] = df["pval"].apply(lambda x: sr.pval_stars(x))
df.loc[:, "N_pos"] = (2 * df["N"] * df["A1FREQ"]).astype(int)
df.head()

# %%
phenos_to_remove = [
    "Severe obesity",
    "Date E66 first reported (obesity)",
    "Date I10 first reported (essential (primary) hypertension)",
    "Myocardial disease incl. angina pectoris",
    "Pulse rate, automated reading",
    "Heart attack diagnosed by doctor",
    "Angina diagnosed by doctor",
    "Fracture resulting from simple fall",
    "Stroke diagnosed by doctor",
    "Cholesterol medication",
]


# %%
plt_df = (
    df.loc[(df.GENE == GENE)]
    .sort_values(by=["Phenotype", "AAF"], ascending=[True, False])
    .groupby(["Phenotype", "MASK"])
    .first()
    .reset_index()
)


plt_df = plt_df.loc[~plt_df.Phenotype.astype(str).isin(phenos_to_remove), :]

effect = {"BT": "OR", "QT": "BETA"}[TRAIT]

group_by_mean = (
    pd.DataFrame({"mean": plt_df.groupby(["Phenotype"]).agg("mean")[effect]})
    .sort_values(by="mean", ascending=False)
    .reset_index()
)

sorter = group_by_mean.Phenotype.tolist()

plt_df.loc[:, "Phenotype"] = plt_df.loc[:, "Phenotype"].astype("category")
plt_df.loc[:, "Phenotypes"] = plt_df.loc[:, "Phenotype"].cat.set_categories(sorter)

plt_df = plt_df.sort_values(
    by=["Phenotype", "MASK"], ascending=[True, False]
).reset_index(drop=True)

phenotypes = plt_df.Phenotype.unique()

print(len(plt_df))
plt_df.head()


# %%
def plot_BT(
    df, width=4, height=12, phenotypes=None, masks=None, lw=1.3, ms=10, fudge=0.6
):

    df = df.loc[df.TRAIT == "BT", :]

    df.loc[:, "CI"] = +(
        "("
        + df.loc[:, "OR_low"].round(2).astype(str)
        + ","
        + df.loc[:, "OR_up"].round(2).astype(str)
        + ")"
    )

    just = df.applymap(lambda x: len(str(x)) + 1).max()

    df.loc[:, "label"] = (
        df.loc[:, "N_pos"].astype(str).str.ljust(just["N_pos"])
        + df.loc[:, "OR"].round(2).astype(str).str.ljust(5)
        + df.loc[:, "CI"].str.ljust(just["CI"])
        + "p= "
        + df.loc[:, "pval"].apply(lambda p: f"{p:.2e}").str.ljust(4)
    )

    if phenotypes is None:
        phenotypes = df.Phenotype.unique()

    num_pheno = len(phenotypes)

    height = 2 * num_pheno

    masks = df.MASK.unique()
    colors = plt.cm.viridis(np.linspace(0, 1, len(masks)))

    fig, axes = plt.subplots(nrows=num_pheno, sharex=True, figsize=(width, height))
    legend_elements = [
        Line2D([0], [0], color=color, marker="o", ms=ms, lw=0, label=mask)
        for color, mask in zip(colors, masks)
    ]

    for i, ax in enumerate(axes):

        temp = df.loc[df.Phenotype.eq(phenotypes[i]), :]
        temp = temp.loc[df.MASK.isin(masks)]

        for mask, color in zip(
            masks,
            colors,
        ):
            temp1 = temp.loc[temp.MASK == mask, :]
            xerr = [temp1["OR_low_lim"].values, temp1["OR_up_lim"].values]

            ax.errorbar(
                temp1["OR"],
                temp1.index,
                alpha=0.99,
                xerr=xerr,
                fmt="o",
                c=color,
                ecolor="black",
                ms=ms,
                mew=0.0,
                mec="black",
                elinewidth=lw,
            )

        ax0 = ax.twinx()
        ax0.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.set_ylim(temp.index[0] - fudge, temp.index[1] + fudge)
        ax.set_yticks([temp.index[0], temp.index[1]])
        ax.set_xlim(0 - fudge, df.OR_up.max() + fudge)
        ax.set_xticks(np.arange(0, int(df.OR_up.max()) + 2, 1))

        ax.set_yticklabels(temp.label, fontsize=9, fontdict={"family": "monospace"})

        ax.tick_params(right=False)
        ax.spines["top"].set_alpha(0)
        ax.spines["left"].set_alpha(0)
        ax.spines["right"].set_alpha(0)
        ax.spines["bottom"].set_alpha(0)

        ax0.tick_params(right=False)
        ax0.tick_params(left=False)
        ax0.spines["top"].set_alpha(0)
        ax0.spines["left"].set_alpha(0)
        ax0.spines["right"].set_alpha(0)

        ax0.grid(False)

        # only show every 3rd yticklabel
        labels = [
            l if i % len(masks) == 0 else "" for i, l in enumerate(temp.Phenotype)
        ]
        ax0.set(yticks=temp.index, yticklabels=labels[::-1])

        ax.axvline(x=1, linestyle="--", color="#4f4f4f")

        ax0.spines["top"].set_alpha(0)
        ax0.spines["left"].set_alpha(0)
        ax0.spines["right"].set_alpha(0)

        if i != len(phenotypes) - 1:
            ax0.spines["bottom"].set_alpha(0)
            ax.tick_params(bottom=False)

        if i == len(phenotypes) - 1:
            ax.set_xlabel("OR")
            ax.tick_params(bottom=True)
            ax.legend(handles=legend_elements, loc="lower right")

    return temp


temp = plot_BT(plt_df)


# %%
def make_plt_df(
    df: pd.DataFrame,
    phenos_to_remove,
    trait=TRAIT,
    gene=GENE,
) -> pd.DataFrame:

    plt_df = (
        df.loc[(df.GENE == GENE)]
        .sort_values(by=["Phenotype", "AAF"], ascending=[True, False])
        .groupby(["Phenotype", "MASK"])
        .first()
        .reset_index()
    )

    plt_df = plt_df.loc[~plt_df.Phenotype.isin(phenos_to_remove), :]

    effect = {"BT": "OR", "QT": "BETA"}[TRAIT]

    group_by_mean = pd.DataFrame(
        {"mean": plt_df.groupby(["Phenotype"]).agg("mean")[effect]}
    ).reset_index()

    group_by_mean = group_by_mean.sort_values(by="mean", ascending=False).reset_index()

    sorter = list(group_by_mean["Phenotype"])

    plt_df.loc[:, "Phenotype"] = plt_df.loc[:, "Phenotype"].astype("category")
    plt_df.loc[:, "Phenotype"].cat.set_categories(sorter, inplace=True)

    plt_df = plt_df.sort_values(
        by=["Phenotype", "MASK"], ascending=[True, False]
    ).reset_index(drop=True)

    phenotypes = plt_df.Phenotype.unique()

    return plt_df


plt_df = make_plt_df(df, phenos_to_remove=phenos_to_remove)
plt_df.head()

# %%
for ax in range(0, len(phenotypes)):
    temp = plt_df.loc[plt_df.Phenotype.eq(phenotypes[ax]), :]
    temp = temp.loc[plt_df.MASK.isin(masks)]

    temps = [temp[temp.MASK == mask] for mask in masks]

    if TRAIT == "BT":
        xerrs = [[temp["OR_se_l"].values, temp["OR_se_u"].values] for temp in temps]

    elif TRAIT == "QT":
        xerrs = [temp["OR_se"] for temp in temps]

    axes[ax].errorbar(
        temps[0]["OR"],
        temps[0].index,
        alpha=0.99,
        xerr=xerrs[0],
        fmt="o",
        c="tab:grey",
        ecolor="black",
        ms=ms,
        mew=0.0,
        mec="black",
        elinewidth=lw,
    )
    axes[ax].errorbar(
        temps[1]["OR"],
        temps[1].index,
        alpha=0.99,
        xerr=xerrs[1],
        fmt="o",
        c="tab:purple",
        ecolor="black",
        ms=ms,
        mew=0.0,
        mec="black",
        elinewidth=lw,
    )
    # axes[ax].errorbar(
    #    temps[2]["OR"],
    #    temps[2].index,
    #    alpha=0.99,
    #    xerr=xerrs[2],
    #    fmt="o",
    #    c="tab:orange",
    #    ecolor="black",
    #    ms=ms,
    #    mew=0.0,
    #    mec="black",
    #    elinewidth=lw,
    # )

    ax0 = axes[ax].twinx()
    ax0.set_ylim([0.25, 3.25])
    ax0.set_yticks(
        [
            1,
            2,
        ]
    )

    y2labels = (
        temp.N_pos.astype(str)
        #    + temp.or_e.astype(str).str.ljust(5).values
        #    + ["("]
        #    + temp.lower_e.values
        #    + [","]
        #    + temp.upper_e.values
        + ["   p = "]
        + temp.pval_e.values
    )

    ax0.set_yticklabels(y2labels, fontsize=9, fontdict={"family": "monospace"})
    ax0.tick_params(right=False)
    ax0.spines["top"].set_alpha(0)
    ax0.spines["left"].set_alpha(0)
    ax0.spines["right"].set_alpha(0)
    ax0.spines["bottom"].set_alpha(0)
    ax0.grid(False)
    # axes[ax].invert_xaxis()
    # only show every 3rd yticklabel
    labels = [l if i % 3 == 0 else "" for i, l in enumerate(temp.Phenotype)]
    axes[ax].set(yticks=temp.index, yticklabels=labels[::-1])

    # axes[ax].axvline(x=0, linestyle="--", color="#4f4f4f")
    axes[ax].tick_params(left=False)

    if TRAIT == "BT":
        axes[ax].set_xlim([-0.2, 5.2])
        axes[ax].axvline(x=1, linestyle=":", color="#4f4f4f")
    elif TRAIT == "QT":
        axes[ax].set_xlim([-0.22, 0.82])
        axes[ax].axvline(x=0, linestyle=":", color="#4f4f4f")

    axes[ax].spines["top"].set_alpha(0)
    axes[ax].spines["left"].set_alpha(0)
    axes[ax].spines["right"].set_alpha(0)
    if ax != len(phenotypes) - 1:
        axes[ax].spines["bottom"].set_alpha(0)
        axes[ax].tick_params(bottom=False)

    if ax == len(phenotypes) - 1:
        axes[ax].set_xlabel(xlab)
        # ax0.legend(handles = legend_elements, loc = "lower right")
#
# axes[len(phenotypes) - 1].set_xticks([0,1,2,3,4,5,6,7,])
for ax in axes.flat:
    ax.margins(0.3)

plt.subplots_adjust(right=1)
# fig.suptitle(title)

# plt.savefig(f"/opt/notebooks/gogoGPCR/tmp/{GENE}_{TRAIT}_new.svg", dpi = 600, bbox_inches="tight", format = "svg")
plt.show()

# %%

# %%
