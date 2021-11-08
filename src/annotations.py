import hail as hl
from src.resources import lauryns_variants
import subprocess

def annotate_OPRM1(mt: hl.MatrixTable, anno_file: str) -> hl.MatrixTable:
    annot = hl.import_table(anno_file, impute = True)

    annot = annot.filter(annot["Transcript stable ID"] == "ENST00000229768")

    annot = annot.annotate(
        protCons=annot["Protein allele"].split("/")[0]
        + hl.str(annot["Variant start in translation (aa)"])
        + annot["Protein allele"].split("/")[1],
        locus = hl.locus("chr" + hl.str(annot["Chromosome/scaffold name"]), annot["Chromosome/scaffold position start (bp)"]),
        alleles = hl.array(annot["Consequence specific allele"].split("/")))
    
    annot = annot.key_by("locus", "alleles")
    
    mt = mt.annotate_rows(annotations=annot[mt.locus, mt.alleles])
    
    lauryns_dict = hl.literal({var: True for var in lauryns_variants})
    
    mt = mt.annotate_rows(is_lauryns = lauryns_dict.get(mt.annotations.protCons))
    
    return mt
    
    
def annotate_GCGR(mt: hl.MatrixTable, anno_file: str) -> hl.MatrixTable:
    ht = hl.import_table(anno_file).key_by("variant")
    mt = mt.annotate_rows(**ht[mt.protCons])
    
    return mt

def annotate_DRD2(mt: hl.MatrixTable, anno_file: str)  -> hl.MatrixTable:
    
    run(["hadoop", "fs", "-put", anno_file, "/tmp/"], check = True, shell = False)
    
    ht = hl.import_table("/tmp/" + anno_file, impute = True).key_by("AA consequence")
    
    mt = mt.annotate_rows(annotations = ht[mt.protCons])

    mt = mt.annotate_rows(Gi1 = mt.annotations.number_of_impairments_Gi1 > 0,
                          GoA = mt.annotations.number_of_impairments_GoA > 0,
                          Gz = mt.annotations.number_of_impairments_Gz > 0)

    mt = mt.annotate_rows(annotations = hl.case()
                          .when(~mt.Gi1 & ~mt.GoA & ~mt.Gz, "WT")
                          .when(mt.Gi1 & ~mt.GoA & ~mt.Gz, "Gi1")
                          .when(~mt.Gi1 & mt.GoA & ~mt.Gz, "GoA")
                          .when(~mt.Gi1 & ~mt.GoA & mt.Gz, "Gz")
                          .when(mt.Gi1 & mt.GoA & ~mt.Gz, "Gi1_GoA")
                          .when(mt.Gi1 & ~mt.GoA & mt.Gz, "Gi1_Gz")
                          .when(~mt.Gi1 & mt.GoA & mt.Gz, "GoA_Gz")
                          .when(mt.Gi1 & mt.GoA & mt.Gz, "Gi1_GoA_Gz")
                         .or_missing())
    return mt