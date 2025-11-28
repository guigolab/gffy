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