# Constants for the GFF parser



SKIP_FEATURES = frozenset(
    {
        "region", 
        "chromosome", 
        "scaffold", 
        "repeat_region", 
        "contig", 
        "scaffold_region", 
        "supercontig",
        "centromere",
        "telomere",
        "biological_region",
    }) 


COMMON_TRANSCRIPT_TYPES = frozenset(
    {
        "mRNA",
        "lncRNA",
        "miRNA",
        "snRNA",
        "snoRNA",
        "scRNA",
        "transcript",
        "rRNA",
        "tRNA",
    }
)

COMMON_GENE_TYPES = frozenset(
    {
        "gene",
        "pseudogene",
        "ncRNA_gene",
    }
)

FEATURE_CODES = {
    "exon": 1,
    "CDS": 2,
}
