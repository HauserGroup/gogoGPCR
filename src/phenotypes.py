import re
import itertools
from warnings import warn


def get_columns_to_keep(df, threshold=200) -> list:
    """
    This function drops all columns which contain more null values than threshold
    :param df: A PySpark DataFrame
    """
    null_counts = (
        df.select([F.count(F.when(~F.col(c).isNull(), c)).alias(c) for c in df.columns])
        .collect()[0]
        .asDict()
    )
    to_keep = [k for k, v in null_counts.items() if v > threshold]
    # df = df.select(*to_keep)

    return to_keep


def new_names(s: str) -> str:
    """
    Fixes a column header for PHESANT use
    """
    s = s.replace("p", "x").replace("i", "")

    if bool(re.search("_\d$", s)):
        s += "_0"
    else:
        s += "_0_0"
    return s


def get_age_sex(participant, fields=["31", "21022"]) -> list:
    """Create a list of age and sex field names

    Parameters
    ----------
    participant : UKB Spark df
        participant df
    fields : list, optional
        list of age and sex fields. Defaults are necessary for PHESANT, by default ["31", "21022"]

    Returns
    -------
    list
        age sex fields
    """
    age_sex = "|".join(fields)
    age_sex_fields = list(
        participant.find_fields(lambda f: bool(re.match(f"^p({age_sex})$", f.name)))
    )

    return [f.name for f in age_sex_fields]


def get_pheno_fields(participant, fields: list) -> list:
    """get a list of phenotype field names for phenotypes of interest

    Parameters
    ----------
    participant : UKB spark df
        participant df
    fields : list
        list of phenotypes of interest, numbered as in showcase

    Returns
    -------
    list
        list of field names
    """
    phenos = "|".join(fields)

    pheno_fields = list(
        participant.find_fields(lambda f: bool(re.match(f"^p({phenos})\D", f.name)))
    )

    return [f.name for f in pheno_fields]


def concatenate(*lists) -> list:
    """concatenate inpout lists to new list

    Returns
    -------
    list
        concatenated list
    """

    return list(itertools.chain(*lists))


def fix_colnames(df):
    """returns df with fixed columns names for PHESANT input

    Parameters
    ----------
    df : UKB Spark df
        from df = participant.retrieve_fields(names=field_names, engine=dxdata.connect())

    Returns
    -------
    Spark df
        with fixed columns
    """

    colnames = ["xeid"] + [new_names(s) for s in df.columns[1:]]

    print(colnames[:10])

    return df.toDF(*colnames)


def filter_to_200k(
    df,
    filter_path="/mnt/project/Data/filters/samples_in_200k_exomes.csv",
    spark=None,
):
    warn("This function is deprecated if using 450k exomes", DeprecationWarning, stacklevel=2)
    """filter df to only samples with WES data

    Parameters
    ----------
    df : UKB spark df
        from participant or similar
    filter_path : str, optional
        location of file with single column of sample ids, for example generated
        by cohort browser, by default "/mnt/project/Data/filters/samples_in_200k_exomes.csv"

    Returns
    -------
    Spark df
        df of length 200k
    """

    samples_with_exomes = spark.read.csv(filter_path, header=True)
    samples_with_exomes = samples_with_exomes.toDF("xeid")

    df = df.join(samples_with_exomes, on="xeid", how="leftsemi")

    print(f"Samples with exomes: {df.count()}")

    return df
