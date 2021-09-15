import pandas as pd
from pathlib import Path
import hail as hl
from typing import List

def get_position(gene: str, mapping: dict):

    blocks = mapping.get("VCF_block")
    blocks = blocks.split(",")

    chromosome = mapping.get("GRCh38_region")
    start = mapping.get("GRCh38_start")
    end = mapping.get("GRCh38_start")

    return chromosome, blocks, start, end

def lookup_vcfs(mapping: dict, vcfdir: str, gene: str, version: str):

    chromosome, blocks, _, _ = get_position(gene, mapping)

    # Locate VCF file(s)

    vcf_files = [
        f"file://{vcfdir}/ukb23156_c{chromosome}_b{block}_{version}.vcf.gz"
        for block in blocks
    ]

    tbi_files = [
        f"file://{vcfdir}/ukb23156_c{chromosome}_b{block}_{version}.vcf.gz.tbi"
        for block in blocks
    ]

    return {"vcfs": vcf_files, "tbis": tbi_files}

def are_variants_present(mt: hl.MatrixTable, var_col: str, variants: List[str]) -> None:
    for v in variants:    
        print(f"{v}: {mt.aggregate_rows(hl.agg.any(mt[var_col] == v))}")