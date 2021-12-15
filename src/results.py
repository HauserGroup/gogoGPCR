from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

def pval_stars(x: float) -> str:
    """[summary]

    Parameters
    ----------
    x : float
        [description]

    Returns
    -------
    str
        [description]
    """
    if x <= 0.0001:
        return "****"
    elif x <= 0.001:
        return "***"
    elif x <= 0.01:
        return "**"
    elif x <= 0.05:
        return "*"
    if x > 0.05:
        return ""


def pheno_search(
    x, ukb_coding: pd.DataFrame, custom_coding: pd.DataFrame
) -> dict:
    """[summary]

    Parameters
    ----------
    x : [type]
        [description]
    ukb_coding : pd.DataFrame
        [description]
    custom_coding : pd.DataFrame
        [description]

    Returns
    -------
    [type]
        [description]
    """
    ukb = dict(zip(ukb_coding.FieldID.astype(str), ukb_coding.Field))
    custom = dict(zip(custom_coding.FieldID.astype(str), custom_coding.Field))
    both = {**ukb, **custom}

    return both[str(x)]

def plot_BT(
    df,
    width=4,
    height=12,
    phenotypes=None,
    masks=None,
    lw=1.3,
    ms=10,
    fudge=0.6,
    xlim=[0, 12],
    title=None,
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
        + "p="
        + df.loc[:, "pval"].apply(lambda p: f"{p:.2e}").str.ljust(4)
        + " "
        + df.loc[:, "pval_stars"].str.ljust(4)
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

    for i, ax in tqdm(enumerate(axes)):

        temp = df.loc[df.Phenotype.eq(phenotypes[i]), :]
        temp = temp.loc[df.MASK.isin(masks)]

        for color, mask in zip(
            colors,
            masks,
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
        print(temp.label)
        ax0 = ax.twinx()
        ax0.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.set_ylim(temp.index[0] - fudge, temp.index[-1] + fudge)
        ax.set_yticks(temp.index)
        # ax.set_xlim(0 - fudge, df.OR_up.max() + fudge)
        ax.set_xlim(xlim[0] - fudge, xlim[1] + fudge)
        # ax.set_xticks(np.arange(0, int(df.OR_up.max()) + 2, 1))
        ax.set_xticks(np.arange(xlim[0], xlim[1] + 1, 1))

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
            ax.legend(handles=legend_elements[::-1], loc="lower right")

    if title:
        fig.suptitle(title)
    plt.subplots_adjust(right=1)

    return fig

def plot_QT(
    df,
    width=4,
    height=None,
    phenotypes=None,
    masks=None,
    lw=1.3,
    ms=10,
    fudge=0.6,
    xlim=[-1, 1],
    title=None,
):

    df = df.loc[df.TRAIT == "QT", :]

    just = df.applymap(lambda x: len(str(x)) + 1).max()

    df.loc[:, "label"] = (
        df.loc[:, "N_pos"].astype(str).str.ljust(just["N_pos"])
        + df.loc[:, "BETA"].round(2).astype(str).str.ljust(5)
        + "± "
        + df.loc[:, "SE"].round(2).astype(str).str.ljust(5)
        + "p="
        + df.loc[:, "pval"].apply(lambda p: f"{p:.2e}").str.ljust(5)
        + " "
        + df.loc[:, "pval_stars"].str.ljust(4)
    )

    if phenotypes is None:
        phenotypes = df.Phenotype.unique()

    num_pheno = len(phenotypes)
    if not height:
        height = 2 * num_pheno

    masks = df.MASK.unique()
    colors = plt.cm.viridis(np.linspace(0, 1, len(masks)))

    fig, axes = plt.subplots(nrows=num_pheno, sharex=True, figsize=(width, height))
    legend_elements = [
        Line2D([0], [0], color=color, marker="o", ms=ms, lw=0, label=mask)
        for color, mask in zip(colors, masks)
    ]

    for i, ax in tqdm(enumerate(axes)):

        temp = df.loc[df.Phenotype.eq(phenotypes[i]), :]
        temp = temp.loc[df.MASK.isin(masks)]

        for color, mask in zip(
            colors,
            masks,
        ):
            temp1 = temp.loc[temp.MASK == mask, :]
            # xerr = [temp1["OR_low_lim"].values, temp1["OR_up_lim"].values]

            ax.errorbar(
                temp1["BETA"],
                temp1.index,
                alpha=0.99,
                xerr=temp1["SE"],
                fmt="o",
                c=color,
                ecolor="black",
                ms=ms,
                mew=0.0,
                mec="black",
                elinewidth=lw,
            )
        print(temp.label)
        ax0 = ax.twinx()
        ax0.yaxis.tick_left()
        ax.yaxis.tick_right()
        ax.set_ylim(temp.index[0] - fudge, temp.index[-1] + fudge)
        ax.set_yticks(temp.index)
        # ax.set_xlim(0 - fudge, df.OR_up.max() + fudge)
        ax.set_xlim(xlim[0] - fudge, xlim[1] + fudge)
        # ax.set_xticks(np.arange(0, int(df.OR_up.max()) + 2, 1))
        ax.set_xticks(np.arange(xlim[0], xlim[1] + 1, 1))

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

        ax.axvline(x=0, linestyle="--", color="#4f4f4f")

        ax0.spines["top"].set_alpha(0)
        ax0.spines["left"].set_alpha(0)
        ax0.spines["right"].set_alpha(0)

        if i != len(phenotypes) - 1:
            ax0.spines["bottom"].set_alpha(0)
            ax.tick_params(bottom=False)

        if i == len(phenotypes) - 1:
            ax.set_xlabel("β")
            ax.tick_params(bottom=True)
            ax.legend(handles=legend_elements[::-1], loc="lower right")

    if title:
        fig.suptitle(title)
    plt.subplots_adjust(right=1)

    return fig