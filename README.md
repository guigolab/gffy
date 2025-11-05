# gffy - GFF3 Genomic Statistics Calculator

[![PyPI version](https://badge.fury.io/py/gffy.svg)](https://badge.fury.io/py/gffy)
[![Python Versions](https://img.shields.io/pypi/pyversions/gffy.svg)](https://pypi.org/project/gffy/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A fast and efficient tool for computing comprehensive statistics from GFF3 genomic annotation files. 

`gffy` processes GFF3 files (including compressed `.gz` files) and generates detailed statistics about genes, transcripts, exons, introns, and CDS features, organized by gene categories (coding genes, non-coding genes, and pseudogenes).

## Features

- 🚀 **High Performance**: Streaming parser that processes large GFF3 files efficiently without loading entire files into memory
- 📊 **Comprehensive Statistics**: Computes detailed metrics for genes, transcripts, exons, introns, and CDS features
- 🧬 **Smart Categorization**: Automatically categorizes genes into coding, non-coding, and pseudogenes
- 🌐 **Remote File Support**: Can directly process GFF3 files from URLs (including `.gz` compressed files)
- 📈 **Detailed Metrics**: Provides count, min, max, mean, and median statistics for various genomic features
- 🔄 **Transcript Type Analysis**: Breaks down statistics by transcript types (mRNA, lnc_RNA, miRNA, etc.)

## Installation

### From PyPI

```bash
pip install gffy
```

### From Source

```bash
git clone https://github.com/emilior/gffy.git
cd gffy
pip install -e .
```

## Quick Start

### Python API

```python
from gffy import compute_gff_stats

# Process a GFF3 file from a URL
url = "https://ftp.ensembl.org/pub/release-110/gff3/homo_sapiens/Homo_sapiens.GRCh38.110.gff3.gz"
stats = compute_gff_stats(url)

# Access statistics
print(f"Coding genes: {stats['coding_genes']['count']}")
print(f"Non-coding genes: {stats['non_coding_genes']['count']}")
print(f"Pseudogenes: {stats['pseudogenes']['count']}")

# Access detailed transcript statistics
for transcript_type, data in stats['coding_genes']['transcripts']['types'].items():
    print(f"{transcript_type}: {data['count']} transcripts")
```

### As a Module

You can also run the included script to process annotations from an API:

```python
from gffy.compute_stats import compute_gff_stats
from gffy.tools.api import get_annotations_without_stats, update_annotation

# Fetch annotations from your API
annotations = get_annotations_without_stats()

# Process each annotation
for annotation in annotations:
    annotation_id = annotation['annotation_id']
    gff_url = annotation['source_file_info']['url_path']
    
    stats = compute_gff_stats(gff_url)
    update_annotation(annotation_id, stats)
```

## Output Structure

The statistics are returned as a JSON-compatible dictionary with the following structure:

```json
{
  "coding_genes": {
    "count": 22178,
    "length_stats": {
      "min": 10,
      "max": 2960899,
      "mean": 48306.64,
      "median": 17025.5
    },
    "transcripts": {
      "count": 102504,
      "per_gene": 4.62,
      "types": {
        "mRNA": {
          "count": 66153,
          "per_gene": 3.05,
          "exons_per_transcript": 9.08,
          "length_stats": { ... },
          "spliced_length_stats": { ... },
          "exon_length_stats": { ... }
        },
        ...
      }
    },
    "features": {
      "exons": { "count": 759800, "length_stats": { ... } },
      "introns": { "count": 657296, "length_stats": { ... } },
      "cds": { "count": 527234, "length_stats": { ... } }
    }
  },
  "non_coding_genes": { ... },
  "pseudogenes": { ... }
}
```

### Gene Categories

- **coding_genes**: Genes with CDS features or protein_coding biotype
- **non_coding_genes**: Genes with exons but no CDS features
- **pseudogenes**: Genes with feature type "pseudogene"

### Statistics Computed

For each category, `gffy` computes:

- **Gene counts and length statistics** (min, max, mean, median)
- **Transcript counts** (total and per-gene average)
- **Per-type transcript statistics**:
  - Count and per-gene ratio
  - Exons per transcript
  - Genomic span length (start to end)
  - Spliced length (sum of exon lengths)
  - Exon length statistics
- **Feature statistics**:
  - Exon counts and lengths
  - Intron counts and lengths
  - CDS counts and lengths (for coding genes)

## Configuration

The package can be configured using environment variables when working with the API integration:

```bash
export API_URL="http://your-api-server:5002"
export AUTH_KEY="your-auth-key"
```

## Use Cases

- **Genome Annotation QC**: Validate and assess the quality of genome annotations
- **Comparative Genomics**: Compare gene structure statistics across different species or assemblies
- **Annotation Pipelines**: Integrate into automated annotation workflows
- **Research**: Analyze transcript diversity, gene structure, and feature distributions

## Performance

`gffy` is optimized for performance:

- Streaming parser (low memory footprint)
- Efficient data structures using `array.array` for coordinates
- String interning for repeated values
- Single-pass processing with orphan resolution

Typical processing time for a complete human genome annotation (GRCh38): ~2-5 minutes on standard hardware.

## Requirements

- Python 3.9+
- requests >= 2.25.0

## JSON Schema

A JSON schema for validating the output is included in `schema.json`. This can be used to validate the statistics output programmatically.

## Documentation

For more detailed information about how statistics are calculated, see [STATS_CALCULATION_GUIDE.md](STATS_CALCULATION_GUIDE.md).

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use `gffy` in your research, please cite:

```
Emilio R. (2025). gffy: Fast GFF3 Genomic Statistics Calculator. 
https://github.com/emilior/gffy
```

## Changelog

### Version 0.1.0 (2025-10-14)

- Initial release
- Support for GFF3 file processing from URLs
- Comprehensive statistics for genes, transcripts, and features
- Categorization by gene type (coding, non-coding, pseudogene)
- Per-transcript-type statistics
- JSON schema for output validation

## Support

For issues, questions, or contributions, please visit the [GitHub repository](https://github.com/emilior/gffy).

