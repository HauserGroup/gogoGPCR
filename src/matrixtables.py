import hail as hl
import annotations as annotations
import utils as utils
import pandas as pd
from functools import partial

from pathlib import Path
from typing import Optional, Union, List


def import_mt(
    genes: List[str],
    mapping: pd.DataFrame,
    vcf_dir: str,
    vcf_version: str,
) -> hl.matrixtable.MatrixTable:
    """Import VCF file or list of VCF files as MatrixTable

    Parameters
    ----------
    vcf_files : Union[str, list[str]]
        VCF or list of VCF files

    Returns
    -------
    hl.matrixtable.MatrixTable
        Raw MatrixTable of all samples and variants, very large. GRCh38 as reference.
    """

    get_vcfs = partial(
        utils.lookup_vcfs, mapping=mapping, vcfdir=vcf_dir, version=vcf_version
    )
    get_regions = partial(utils.lookup_regions, mapping=mapping)

    # evil double list-comprehension
    vcf_files = [vcf for gene in genes for vcf in get_vcfs(gene=gene)]
    regions = [region for gene in genes for region in get_regions(gene=gene)]

    mts = hl.import_gvcfs(
        vcf_files,
        partitions=regions,
        reference_genome="GRCh38",
        array_elements_required=False,
    )

    if len(mts) == 1:
        return mts[0]
    else:
        return hl.MatrixTable.union_rows(*mts)


def downsample_mt(
    mt: hl.matrixtable.MatrixTable,
    prob: Optional[float] = None,
    seed=42,
) -> hl.matrixtable.MatrixTable:
    """Reduce number of samples in MatrixTable

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        MatrixTable with samples as columns
    prob : Optional[float], optional
        Fraction of samples to keep, 1 / 200 results in ~1000 samples, by default None

    Returns
    -------
    hl.matrixtable.MatrixTable
        Downsampled MatrixTable
    """

    if prob is not None:
        return mt.sample_cols(p=prob, seed=seed)
    else:
        return mt


def interval_qc_mt(
    mt: hl.matrixtable.MatrixTable,
    bed_file: Union[str, Path],
) -> hl.matrixtable.MatrixTable:
    """Filter to only Target region used by the WES capture experiment

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        MatrixTable
    intervals : str
        .BED file of targeted capture regions which meet quality standards

    Returns
    -------
    hl.matrixtable.MatrixTable
        MatrixTable filtered to only target regions
    """

    interval_table = hl.import_bed(
        bed_file,
        reference_genome="GRCh38",
    )

    mt = mt.filter_rows(hl.is_defined(interval_table[mt.locus]))

    return mt


def sample_QC_mt(
    mt: hl.matrixtable.MatrixTable,
    MIN_CALL_RATE: float,
    MIN_MEAN_DP: float,
    MIN_MEAN_GQ: float,
) -> hl.matrixtable.MatrixTable:
    """Filter out samples failing quality control

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        Matrix table with sample information
    samples_to_remove : str
        List of samples to filter out as generated by sample_hard_filter.py
        includes relatedness, sex aneuploidy, outliers etc.

    Returns
    -------
    hl.matrixtable.MatrixTable
        MatrixTable with filtered samples
    """
    mt = hl.sample_qc(mt)
    mt = mt.filter_cols(
        (mt.sample_qc.call_rate >= MIN_CALL_RATE)
        & (mt.sample_qc.dp_stats.mean >= MIN_MEAN_DP)
        & (mt.sample_qc.gq_stats.mean >= MIN_MEAN_GQ)
    )

    return mt


def smart_split_multi_mt(
    mt: hl.matrixtable.MatrixTable, left_aligned=False
) -> hl.matrixtable.MatrixTable:

    mt = mt.key_rows_by("locus", "alleles")

    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    split = hl.split_multi_hts(multi, left_aligned=left_aligned)
    mt = split.union_rows(bi)

    return mt


def variant_QC_mt(
    mt: hl.matrixtable.MatrixTable,
    MIN_P_HWE: Optional[float],
    MIN_GQ: Union[float, int, None],
) -> hl.matrixtable.MatrixTable:

    mt = hl.variant_qc(mt)

    if MIN_P_HWE is not None:
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe >= MIN_P_HWE)

    if MIN_GQ is not None:
        mt = mt.filter_rows(mt.variant_qc.gq_stats.mean >= MIN_GQ)

    mt = mt.filter_rows(
        (mt.variant_qc.AF[0] > 0.0) & (mt.variant_qc.AF[0] < 1.0)
    )

    return mt


def genotype_filter_mt(
    mt: hl.matrixtable.MatrixTable,
    MIN_DP: Union[float, int, None],
    MIN_GQ: Union[float, int, None],
    log_entries_filtered: True,
) -> hl.matrixtable.MatrixTable:

    mt = mt.annotate_entries(AB=(mt.AD[1] / hl.sum(mt.AD)))

    # set filter condition for AB
    filter_condition_ab = (
        (mt.GT.is_hom_ref() & (mt.AB <= 0.1))
        | (mt.GT.is_het() & (mt.AB >= 0.25) & (mt.AB <= 0.75))
        | (mt.GT.is_hom_var() & (mt.AB >= 0.9))
    )
    fraction_filtered = mt.aggregate_entries(
        hl.agg.fraction(~filter_condition_ab)
    )

    print(
        f"Filtering {fraction_filtered * 100:.2f}% entries out of downstream"
        " analysis."
    )

    mt = mt.filter_entries(
        (mt.GQ >= MIN_GQ)
        & (mt.DP >= MIN_DP)
        & (
            (mt.GT.is_hom_ref() & (mt.AB <= 0.1))
            | (mt.GT.is_het() & (mt.AB >= 0.25) & (mt.AB <= 0.75))
            | (mt.GT.is_hom_var() & (mt.AB >= 0.9))
        )
    )

    if log_entries_filtered:
        mt = mt.compute_entry_filter_stats()

    return mt


def add_varid(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Annotate rows with varid

    Parameters
    ----------
    mt : hl.MatrixTable
        [description]

    Returns
    -------
    hl.MatrixTable
        [description]
    """

    mt = mt.annotate_rows(
        varid=hl.delimit(
            [
                mt.locus.contig,
                hl.str(mt.locus.position),
                mt.alleles[0],
                mt.alleles[1],
            ],
            ":",
        )
    )

    return mt


#def annotate_mt(mt: hl.MatrixTable, name: str, **kwargs) -> hl.MatrixTable:
#
#    func = getattr(annotations, f"annotate_{name}")
#
#    return func(mt, **kwargs)


def filter_related_mt(
    mt: hl.MatrixTable, rel_file: str, max_kinship: Optional[float]
) -> hl.MatrixTable:
    rel = hl.import_table(
        rel_file,
        delimiter=" ",
        impute=True,
        types={"ID1": "str", "ID2": "str"},
    )

    rel = rel.filter(rel.Kinship > max_kinship, keep=True)

    related_samples_to_remove = hl.maximal_independent_set(
        i=rel.ID1,
        j=rel.ID2,
        keep=False,
    ).key_by("node")

    return mt.anti_join_cols(related_samples_to_remove)


def recode_GT_to_GP(
    mt: hl.matrixtable.MatrixTable,
) -> hl.matrixtable.MatrixTable:

    GPs = hl.literal([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    mt = mt.annotate_entries(GP=GPs[mt.GT.n_alt_alleles()])

    return mt


def write_bgen(mt: hl.matrixtable.MatrixTable, output: str) -> None:

    mt = add_varid(mt)

    mt = recode_GT_to_GP(mt)

    hl.export_bgen(
        mt=mt, varid=mt.varid, rsid=mt.varid, gp=mt.GP, output=output
    )


def generate_report(mt):
    intr = mt.filter_rows((hl.is_defined(mt.annotations)))
    intr = hl.variant_qc(intr)
    intr = intr.select_rows(
        intr.variant_qc, intr.protCons, intr.annotations
    ).rows()
    intr = intr.annotate(**intr.variant_qc)
    intr = intr.annotate(**intr.annotations)
    intr = intr.drop("variant_qc", "gq_stats", "dp_stats", "annotations")

    return intr
