import pandas as pd
from pathlib import Path
import hail as hl
from typing import List, Tuple
from distutils.version import LooseVersion


def get_position(gene: str, mapping: pd.DataFrame):
    """[summary]

    Parameters
    ----------
    gene : str
        [description]
    mapping : pd.DataFrame
        [description]

    Returns
    -------
    [type]
        [description]
    """
    blocks = mapping.loc[gene, "VCF_block"].split(",")
    chromosome = mapping.loc[gene, "GRCh38_region"]
    start = mapping.loc[gene, "GRCh38_start"]
    end = mapping.loc[gene, "GRCh38_end"]

    return chromosome, blocks, start, end


def lookup_regions(gene: str, mapping: pd.DataFrame) -> hl.expr.LocusExpression:
    """[summary]

    Parameters
    ----------
    gene : str
        [description]
    mapping : pd.DataFrame
        [description]

    Returns
    -------
    hl.expr.LocusExpression
        [description]
    """
    chromosome, _, start, end = get_position(gene, mapping)

    region = [
        hl.parse_locus_interval(
            f"[chr{chromosome}:{start}-chr{chromosome}:{end}]"
        )
    ]

    return region


def lookup_vcfs(
    mapping: pd.DataFrame, vcfdir: str, gene: str, version: str, field_id: int = 23148
) -> List[str]:
    """[summary]

    Parameters
    ----------
    mapping : pd.DataFrame
        [description]
    vcfdir : str
        [description]
    gene : str
        [description]
    version : str
        [description]

    Returns
    -------
    List[str]
        [description]
    """

    chromosome, blocks, _, _ = get_position(gene, mapping)

    vcf_files = [
        f"file://{vcfdir}/ukb{field_id}_c{chromosome}_b{block}_{version}.vcf.gz"
        for block in blocks
    ]

    return vcf_files


def are_variants_present(
    mt: hl.MatrixTable, var_col: str, variants: List[str]
) -> None:
    """[summary]

    Parameters
    ----------
    mt : hl.MatrixTable
        [description]
    var_col : str
        [description]
    variants : List[str]
        [description]
    """
    for v in variants:
        print(f"{v}: {mt.aggregate_rows(hl.agg.any(mt[var_col] == v))}")


def get_stats(
    mt: hl.matrixtable.MatrixTable,
) -> Tuple[hl.table.Table, hl.table.Table]:
    """[summary]

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        [description]

    Returns
    -------
    Tuple[hl.table.Table, hl.table.Table]
        [description]
    """
    intr = mt.filter_rows((hl.is_defined(mt.labels)))
    intr = hl.variant_qc(intr)
    intr = intr.rows()
    intr = intr.select(intr.variant_qc, intr.protCons, intr.labels)
    intr = intr.annotate(**intr.variant_qc)
    intr = intr.drop(
        "variant_qc",
        "gq_stats",
        "dp_stats",
    )
    stats = intr.group_by(intr.labels).aggregate(
        n_carriers=hl.agg.sum(intr.n_het), n_variants=hl.agg.count()
    )
    return stats, intr


def haplotype_carriers(mt: hl.MatrixTable, variants: List[str]) -> int:
    """[summary]

    Parameters
    ----------
    mt : hl.MatrixTable
        [description]
    variants : List[str]
        [description]

    Returns
    -------
    int
        [description]
    """
    mt = mt.filter_rows(hl.literal(variants).contains(mt.protCons))

    assert mt.count_rows() == len(
        variants
    ), "Haplotype does not exist, variant missing"

    mt = mt.annotate_cols(carrier=hl.agg.all(mt.GT.is_non_ref()))
    num_carriers = mt.aggregate_cols(hl.agg.sum(mt.carrier))

    return num_carriers


def fields_for_id(field_id, participant):
    """[summary]

    Parameters
    ----------
    field_id : [type]
        [description]
    participant : [type]
        [description]

    Returns
    -------
    [type]
        [description]
    """

    field_id = str(field_id)
    fields = participant.find_fields(
        name_regex=r"^p{}(_i\d+)?(_a\d+)?$".format(field_id)
    )

    return sorted(fields, key=lambda f: LooseVersion(f.name))
