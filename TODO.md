
# **📌 GFF Statistics Parser – TODO List**

## **1. Core Statistics Functions**

* [ ] `total_records()` – Return total number of features
* [ ] `feature_counts()` – Counts of each feature type (gene, exon, CDS, etc.)
* [ ] `biotype_counts()` – Counts of features grouped by biotype
* [ ] `mean_length(feature=None, biotype=None)` – Average feature length (optionally filtered)
* [ ] `median_length(feature=None, biotype=None)` – Median feature length
* [ ] `min_length(feature=None, biotype=None)` – Minimum feature length
* [ ] `max_length(feature=None, biotype=None)` – Maximum feature length
* [ ] `spliced_length(transcript_id=None)` – Sum of exon lengths per transcript
* [ ] `length_distribution(feature=None, biotype=None)` – Return list of lengths for plotting or analysis

---

## **2. Custom Filtering / Query Functions**

* [ ] `filter_records(feature=None, biotype=None, seqid=None)` – Return filtered records
* [ ] `count(feature=None, biotype=None, seqid=None)` – Count filtered records
* [ ] `mean_filtered_length(feature=None, biotype=None, seqid=None)` – Mean length of filtered records
* [ ] `median_filtered_length(feature=None, biotype=None, seqid=None)` – Median length of filtered records
* [ ] `spliced_length_filtered(feature="transcript", biotype=None)` – Spliced length for filtered transcripts
* [ ] `grouped_counts(attributes: list[str])` – Counts grouped by attribute combinations
* [ ] `grouped_mean_length(attributes: list[str])` – Mean lengths grouped by attribute combinations

---

## **3. Quality Control / Sanity Checks**

* [ ] `check_zero_length_features()` – Identify features with zero length
* [ ] `check_single_exon_genes()` – Identify mono-exonic genes
* [ ] `check_missing_cds()` – Identify transcripts missing CDS
* [ ] `check_invalid_coordinates()` – Check start > end or invalid strands
* [ ] `summary_report()` – Generate combined QC report (counts, lengths, suspicious features)

---

## **4. Comparative / Benchmarking Functions**

* [ ] `compare_feature_counts(other_parser)` – Compare counts per feature type between two annotations
* [ ] `compare_length_statistics(other_parser)` – Compare mean/median/min/max lengths
* [ ] `compare_spliced_lengths(other_parser)` – Compare transcript spliced lengths
* [ ] `generate_comparison_table(other_parser, attributes=["feature", "biotype"])` – Combined comparison table

---

## **5. Convenience / Utility Functions**

* [ ] `to_dataframe()` – Export records as a pandas DataFrame
* [ ] `histogram(feature=None, biotype=None, bins=50)` – Return histogram data of lengths
* [ ] `plot_length_distribution(feature=None, biotype=None)` – Plot histogram of lengths
* [ ] `export_summary(path, format="csv")` – Save summary statistics to file

---

## **6. Input Handling / Streaming Utilities**

* [ ] `_open_stream(source)` – Internal: open file-like object (path, URL, gzip)
* [ ] `from_url(url)` – Classmethod: instantiate parser from URL transparently
* [ ] `from_file(path)` – Classmethod: instantiate parser from local file
* [ ] `__enter__` / `__exit__` – Context manager: download URL to temp file if needed, cleanup on exit

---

## **7. Internal / Parsing Helpers**

* [ ] `_parse_stream()` – Iterate over file lines and build records
* [ ] `_parse_line(line)` – Parse a single GFF line into `GFFRecord`
* [ ] `_is_url(s)` – Utility to check if string is URL
* [ ] `_prepare_file(source)` – Download URL to temp file or return local path
* [ ] `_compute_spliced_length()` – Helper to calculate spliced transcript lengths
* [ ] `_cached_stats` – Optional caching mechanism for computed statistics

---

This list **groups all functionality logically**, so you can implement it module by module, starting with **core statistics** and **streaming parsing**, then moving on to **filters, QC, comparisons, and convenience functions**.

---

If you want, I can **turn this TODO list into a ready-to-implement Python class skeleton** with **all methods stubbed and docstrings included**, so you can start coding immediately.

Do you want me to do that?
