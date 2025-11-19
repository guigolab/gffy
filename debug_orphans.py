#!/usr/bin/env python3
"""Debug script to understand orphan resolution."""

from src.gffy.compute_stats import compute_gff_stats

gff_url = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/GCA_000001405.29/ensembl/geneset/2024_11/genes.gff3.gz"

# Monkey patch to add debugging
original_process_feature = None

def debug_compute():
    from src.gffy.compute_stats import process_feature
    from src.gffy.tools.helpers import process_feature as original
    
    global original_process_feature
    original_process_feature = original
    
    call_count = [0]
    processed_count = [0]
    orphan_count = [0]
    
    def debug_wrapper(*args, **kwargs):
        call_count[0] += 1
        result = original(*args, **kwargs)
        if result:
            processed_count[0] += 1
        else:
            if args[6]:  # parent_ids
                orphan_count[0] += 1
        if call_count[0] % 100000 == 0:
            print(f"Processed {call_count[0]} features: {processed_count[0]} successful, {orphan_count[0]} orphans")
        return result
    
    # This won't work easily, let's just run it and check the results
    print("Running compute_gff_stats...")
    result = compute_gff_stats(gff_url)
    
    print("\nChecking results...")
    if result:
        total_genes = sum(stats.get('count', 0) for stats in result.values() if isinstance(stats, dict))
        print(f"Total genes processed: {total_genes}")

if __name__ == "__main__":
    debug_compute()

