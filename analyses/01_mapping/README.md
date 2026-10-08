# 01. Mapping

Split into two pipelines so raw sequence data can be fetched and mapped in disk-bounded
batches instead of all at once — see [`code/01_mapping/README.md`](../../code/01_mapping/README.md)
for why, and [`../00_fetch-raw-data/run_batched.sh`](../00_fetch-raw-data/run_batched.sh) for the
driver that ties them together.

| Pipeline | Directory | Runs |
|---|---|---|
| Fetch + map | [`1_fetch_map/`](1_fetch_map/) | Once per batch (batchable — see `run_batched.sh`) |
| Cohort genotyping | [`2_cohort_genotyping/`](2_cohort_genotyping/) | Exactly once, after every batch has completed |

## Running

Don't run either pipeline directly for a full cohort — use the batch driver:

```bash
export NXF_PROFILE=uppmax
../00_fetch-raw-data/run_batched.sh
```

`run_batched.sh` invokes `1_fetch_map/run_nextflow.sh` once per batch (materializing each
batch's BAM/gVCF/depth files before cleaning up that batch's `work/` directory) and then
`2_cohort_genotyping/run_nextflow.sh` once, over the complete accumulated cohort.

Each pipeline's own `run_nextflow.sh` can still be run directly for a single batch, outside the
batch driver — see each subdirectory's own inputs. **No working stub test currently exists for
either substage** (verified 2026-09-18 — see `code/01_mapping/README.md`'s "Known characteristics"
for what's missing and why); real-data validation happens on the cluster instead.
