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
    """Maps a (list of) genes to their corresponding VCF file and region before import as Hail MatrixTable

    Parameters
    ----------
    genes : List[str]
        List of HGNC genes
    mapping : pd.DataFrame
        DataFrame with columns: "HGNC", "GRCh38_end", "GRCh38_start", "GRCh38_region", "VCF_block". See data folder
    vcf_dir : str
        Location of UKB VCF files
    vcf_version : str
        Used to construct file names, currently v1

    Returns
    -------
    hl.matrixtable.MatrixTable
        MT of all samples with all variants located within specified genes
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

    # union multiple matrixtables before returning them

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
    MIN_CALL_RATE: 0.97,
    MIN_MEAN_DP: 10,
    MIN_MEAN_GQ: 50,
) -> hl.matrixtable.MatrixTable:
    """Add Sample QC and filter out samples failing quality control

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        input MT
    MIN_CALL_RATE : float
        Minimum call rate, default is pretty standard
    MIN_MEAN_DP : float
        Minimum call rate, default is a little lax for exomes
    MIN_MEAN_GQ : float
        Minimum call rate, default is reasonable

    Returns
    -------
    hl.matrixtable.MatrixTable
        MT with sample QC columns fields (.sample_qc) and low-quality samples filtered out
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
    """Split multiple alleles into bi-allelic in a clever way

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        MT with non-bi-allelic sites
    left_aligned : bool, optional
        Assume that alleles are left-aligned for faster splitting, by default False

    Returns
    -------
    hl.matrixtable.MatrixTable
        MT with only bi-allelic sites
    """

    mt = mt.key_rows_by("locus", "alleles")

    # Only split relevant alleles as suggested by Hail docs

    bi = mt.filter_rows(hl.len(mt.alleles) == 2)
    bi = bi.annotate_rows(a_index=1, was_split=False)
    multi = mt.filter_rows(hl.len(mt.alleles) > 2)
    split = hl.split_multi_hts(multi, left_aligned=left_aligned)
    mt = split.union_rows(bi)

    return mt


def variant_QC_mt(
    mt: hl.matrixtable.MatrixTable,
    MIN_P_HWE: Optional[float] = 10 ** -15,
    MIN_GQ: Union[float, int, None] = 20,
) -> hl.matrixtable.MatrixTable:
    """Add variant QC and filter out variants failing (very basic) quality control

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        Input MT
    MIN_P_HWE : Optional[float], optional
        If things went wrong, this is proabbly a reasonable setting, by default 10**-15
    MIN_GQ : Union[float, int, None], optional
        Minimum GQ, by default 20

    Returns
    -------
    hl.matrixtable.MatrixTable
        MT with variant QC row fields (.variant_qc) and without low-quality variants
    """

    mt = hl.variant_qc(mt)

    if MIN_P_HWE is not None:
        mt = mt.filter_rows(mt.variant_qc.p_value_hwe >= MIN_P_HWE)

    if MIN_GQ is not None:
        mt = mt.filter_rows(mt.variant_qc.gq_stats.mean >= MIN_GQ)

    mt = mt.filter_rows((mt.variant_qc.AF[0] > 0.0) & (mt.variant_qc.AF[0] < 1.0))

    return mt


def genotype_filter_mt(
    mt: hl.matrixtable.MatrixTable,
    MIN_DP: Union[float, int, None] = 10,
    MIN_GQ: Union[float, int, None] = 20,
    log_entries_filtered=True,
) -> hl.matrixtable.MatrixTable:
    """[summary]

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        Input MT
    MIN_DP : Union[float, int, None], optional
        Minimum depth, lax default, by default 10
    MIN_GQ : Union[float, int, None], optional
        Minimum quality, lax default, by default 20
    log_entries_filtered : bool, optional
        whether to run compute_entry_filter_stats(), by default True

    Returns
    -------
    hl.matrixtable.MatrixTable
        MT with genotype QC
    """
    mt = mt.annotate_entries(AB=(mt.AD[1] / hl.sum(mt.AD)))

    # set filter condition for AB

    filter_condition_ab = (
        (mt.GT.is_hom_ref() & (mt.AB <= 0.1))
        | (mt.GT.is_het() & (mt.AB >= 0.25) & (mt.AB <= 0.75))
        | (mt.GT.is_hom_var() & (mt.AB >= 0.9))
    )
    fraction_filtered = mt.aggregate_entries(hl.agg.fraction(~filter_condition_ab))

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
    """Annotate rows with varid for export_bgen

    Parameters
    ----------
    mt : hl.MatrixTable
        MT

    Returns
    -------
    hl.MatrixTable
        MT with varid row field
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


# DEPRECATED
# def annotate_mt(mt: hl.MatrixTable, name: str, **kwargs) -> hl.MatrixTable:
#
#    func = getattr(annotations, f"annotate_{name}")
#
#    return func(mt, **kwargs)


def recode_GT_to_GP(
    mt: hl.matrixtable.MatrixTable,
) -> hl.matrixtable.MatrixTable:
    """Recode GT to fake GPs for export_bgen

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        Input MT

    Returns
    -------
    hl.matrixtable.MatrixTable
        MT with GP entry fields
    """

    GPs = hl.literal([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

    mt = mt.annotate_entries(GP=GPs[mt.GT.n_alt_alleles()])

    return mt


def write_bgen(mt: hl.matrixtable.MatrixTable, output: str) -> None:
    """Export MatrixTable as .bgen and .sample file

    Parameters
    ----------
    mt : hl.matrixtable.MatrixTable
        Input MT
    output : str
        Name for output.bgen and output.sample
    """
    mt = add_varid(mt)

    mt = recode_GT_to_GP(mt)

    hl.export_bgen(mt=mt, varid=mt.varid, rsid=mt.varid, gp=mt.GP, output=output)
