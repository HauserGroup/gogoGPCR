import re
import itertools

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
    age_sex = "|".join(fields)
    age_sex_fields = list(
        participant.find_fields(lambda f: bool(re.match(f"^p({age_sex})$", f.name)))
    )
    
    return [f.name for f in age_sex_fields]

def get_pheno_fields(participant, fields) -> list:
    phenos = "|".join(fields)
    
    pheno_fields = list(
        participant.find_fields(lambda f: bool(re.match(f"^p({phenos})\D", f.name)))
    )
    
    return [f.name for f in pheno_fields]

def concatenate(*lists):
    
    return list(itertools.chain(*lists))

def fix_colnames(df):
    
    colnames = ["xeid"] + [new_names(s) for s in df.columns[1:]]

    print(colnames[:10])
    
    return df.toDF(*colnames)

def filter_to_200k(df, filter_path = "/mnt/project/Data/filters/samples_in_200k_exomes.csv"):
    
    samples_with_exomes = spark.read.csv(filter_path, header = True)
    samples_with_exomes = samples_with_exomes.toDF("xeid")
    
    df = df.join(samples_with_exomes, on = "xeid", how = "leftsemi")
    
    print(f"Samples with exomes: {df.count()}")
    
    return df