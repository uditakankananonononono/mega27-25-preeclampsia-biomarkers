# Source provenance and rights documentation gaps (preeclampsia child, 2026-10-08 audit record)

This file records what this repository's own documents say about source hashes and terms. It is a documentation record, not a license verdict, not an integrity certificate, and not a clearance of any recorded rights hold. A matching claim below means two recorded hash strings are equal at a named field in pinned repository documents. It does not mean any assay bytes were re-downloaded or verified. No rights record was found for these sources, so reuse and redistribution are NOT cleared by this file.

Audited child commit: `f938f61fde6aa70ca8e160718444b67314d07122`  
Compared shared-parent repo: `mega27-25-biomarkers-underserved-diseases` at `36ff87995f4fe0e1d08f9e2b2cd885d341b88937`

Rights: No source-specific rights check found in inspected disease-child documentation. Hash matches are not legal clearance or byte verification.

## Per-child correction

Correction recorded by the mapping audit: the hashes in `results/pe_p21_raw_checksums.csv` (GSE148241, 39 rows) are for deposited per-sample raw count files (`.raw.counts.txt.gz`), not GSM pages.

## 1. Result-file hash records in this child (12)

Hash values are shown as 12-hex display prefixes only; the full value is in the cited file and field. Classes are kept separate on purpose: an expression/assay source hash, a GEO series-matrix hash (metadata that may include expression), a reference annotation, and a published artifact are different kinds of record.

- `results/pe_p24_result.json` `/matrix_sha256` hash `db778a7e3094...` (GSE295760); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_p24_result.json` `/workbook_sha256` hash `bb9bfde3c93a...` (GSE295760); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_p14_hdp_reuse_diagnostic.json` `/count_sha256` hash `d9ba5a9d6832...` (GSE186257,GSE234729,GSE303463); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_p23_result.json` `/matrix_sha256` hash `f171107efc41...` (GSE44711); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_p42_result.json` `/source_sha256` hash `f3968cba1cb9...` (GSE284329); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_p44_result.json` `/source_sha256` hash `05fc9d8d9235...` (GSE279757); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_p22_result.json` `/workbook_sha256` hash `c882ad47d8f6...` (GSE177049); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_p19_result.json` `/workbook_sha256` hash `42e900944734...` (GSE143953); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_plac_p13_GSE114691.json` `/source_sha256/ControlONLY` hash `3d429c4f8c43...` (GSE114691); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_plac_p13_GSE114691.json` `/source_sha256/PEONLY` hash `db15986421b4...` (GSE114691); class: expression_or_assay_source_hash_record; equal to shared-parent field
- `results/pe_plac_p12.json` `/sha256` hash `d2c794f403c5...` (GSE190971); class: source_hash_record_target_needs_context; equal to shared-parent field
- `results/pe_plac_p13_GSE186257.json` `/sha256` hash `5ec95f4efab2...` (GSE186257); class: source_hash_record_target_needs_context; equal to shared-parent field

## 2. Additional per-file hash table

- `results/pe_p21_raw_checksums.csv` (GSE148241): 39 rows; class: 39_per_GSM_raw_count_gzip_files; rows equal to shared-parent rows: True. Shared scripts/pe_p21_test.py finds per-GSM .raw.counts.txt.gz URLs and downloads/hash-checks selected39 files; not GSM-page hashes.

## Coverage boundary

- Source accessions listed in the child manifest: 34.
- Of those, 26 have no mapped result-file payload hash in this audit: GSE10588, GSE147776, GSE149440, GSE186819, GSE190639, GSE192902, GSE204835, GSE24129, GSE25906, GSE272342, GSE296973, GSE30186, GSE303840, GSE306864, GSE35574, GSE43942, GSE48424, GSE54400, GSE54618, GSE60438, GSE66273, GSE73374, GSE75010, GSE91189, GSE93839, GSE94643.
- 0 hash fields inherited from other diseases' records are not counted toward this child.
- Records absent from child manifest can still have result-file hashes. Neither presence nor absence proves assay acquisition by this scout. Other-disease copied hashes never count toward this child.
