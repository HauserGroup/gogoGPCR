import pandas as pd


def pval_stars(x: float):
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


def pheno_search(x, ukb_coding: pd.DataFrame, custom_coding: pd.DataFrame):
    ukb = dict(zip(ukb_coding.FieldID.astype(str), ukb_coding.Field))
    custom = dict(zip(custom_coding.FieldID.astype(str), custom_coding.Field))
    both = {**ukb, **custom}

    return both[str(x)]
