#!/usr/bin/env python3
"""Benchmark script to measure runtime and memory usage of compute_gff_stats."""

import tracemalloc
import time
import sys
from src.gffy.compute_stats import compute_gff_stats
import json

def format_bytes(bytes_val: int) -> str:
    """Format bytes to human-readable format."""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.2f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.2f} TB"

def main():
    gff_url = "https://ftp.ebi.ac.uk/pub/ensemblorganisms/Homo_sapiens/GCA_000001405.29/ensembl/geneset/2024_11/genes.gff3.gz"
    
    print("=" * 80)
    print("Benchmarking compute_gff_stats")
    print("=" * 80)
    print(f"GFF Source: {gff_url}")
    print()
    
    # Start memory tracing
    tracemalloc.start()
    
    # Start time measurement
    start_time = time.time()
    
    try:
        # Run the function
        print("Starting computation...")
        result = compute_gff_stats(gff_url)
        
        # End time measurement
        end_time = time.time()
        elapsed_time = end_time - start_time
        
        # Get memory statistics
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Print results
        print()
        print("=" * 80)
        print("RESULTS")
        print("=" * 80)
        print(f"Runtime: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        print(f"Peak Memory Usage: {format_bytes(peak)}")
        print(f"Current Memory Usage: {format_bytes(current)}")
        print()
        
        # Print summary of results
        if result:
            print("Statistics Summary:")
            for category, stats in result.items():
                if stats and "count" in stats:
                    print(f"  {category}: {stats.get('count', 0)} genes")
            with open("stats_refactored.json", "w") as f:
                json.dump(result, f, indent=4)
        else:
            print("No results returned (empty dictionary)")
            
    except Exception as e:
        tracemalloc.stop()
        print(f"\nERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()

