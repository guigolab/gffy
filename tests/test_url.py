import time
import tracemalloc

from gffy import compute_gff_stats


def test_url():
    url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gff.gz"
    
    # Start tracking memory
    tracemalloc.start()
    start_time = time.perf_counter()
    
    # Run the function
    compute_gff_stats(url)
    
    # Measure elapsed time
    end_time = time.perf_counter()
    runtime = end_time - start_time
    
    # Get memory statistics
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    
    # Print measurements
    print(f"\n{'='*60}")
    print(f"Performance Metrics for compute_gff_stats")
    print(f"{'='*60}")
    print(f"Runtime: {runtime:.2f} seconds ({runtime/60:.2f} minutes)")
    print(f"Current memory usage: {current / 1024 / 1024:.2f} MB")
    print(f"Peak memory usage: {peak / 1024 / 1024:.2f} MB")
    print(f"{'='*60}\n")
    
if __name__ == "__main__":
    test_url()