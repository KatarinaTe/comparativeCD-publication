# data/01_mapping

Sample lists and small reference files for [`code/01_mapping`](../../code/01_mapping). The
original pre-conversion `README_01.Mapping.md` (describing the old shell-script-based workflow)
has been moved to
[`archive/legacy_readmes/01_mapping_README_01.Mapping.md`](../../archive/legacy_readmes/01_mapping_README_01.Mapping.md) —
superseded by `code/01_mapping/README.md`.

| File | Description |
|---|---|
| `High_pass_samples.txt` | 11 HighPass sample IDs |
| `All.DA_LowPass.Samples_fastq.txt` | 3,044 LowPass sample IDs (pre-`run_batched.sh` full list) |
| `All.DA_LowPass.Samples_bam.txt` | 3,044 LowPass sample IDs with BAM-file links |
| `All.DA_LowPass.Samples_MoveList` | 3,044 LowPass sample IDs |
| `CHR1-38.txt` | List of the 38 dog autosomes |
| `cohort.sample_map` | HighPass sampleID→gVCF-path map, input to joint SNP calling |
| `exclude_12_LowPass_samples.txt` | 12 samples included through mapping but excluded downstream (no Darwin's Ark dog ID/phenotypes, not in SRA) — consumed by `02_data_processing/3_Merging_filtering`'s `APPLY_EXCLUDE_LOWPASS12` |
