import pandas as pd
from pathlib import Path

def get_position(gene: str, mapping: dict):

    blocks = mapping.get("VCF_block")
    blocks = blocks.split(",")

    chromosome = mapping.get("GRCh38_region")
    start = mapping.get("GRCh38_start")
    end = mapping.get("GRCh38_start")

    return chromosome, blocks, start, end

def lookup_vcfs(mapping: dict, vcfdir: Path, gene: str, version: str):

    chromosome, blocks, _, _ = get_position(gene, mapping)

    # Locate VCF file(s)

    vcf_files = [
        f"{vcfdir}/ukb23156_c{chromosome}_b{block}_{version}.vcf.bgz"
        for block in blocks
    ]

    tbi_files = [
        f"{vcfdir}/ukb23156_c{chromosome}_b{block}_{version}.vcf.bgz.tbi"
        for block in blocks
    ]

    return {"vcfs": vcf_files, "tbis": tbi_files}