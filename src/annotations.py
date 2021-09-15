import hail as hl
from src.resources import lauryns_variants

def annotate_OPRM1(mt: hl.MatrixTable) -> hl.MatrixTable:
    annot = hl.import_table("file:///mnt/project/misc_data/martquery_0901114913_896.txt", impute = True)

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
    
    