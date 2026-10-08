# Container Derivation for 01 Mapping Processes

This directory contains Conda environment specifications used to freeze Docker containers for the Nextflow processes defined in [`code/01_mapping/`](../../code/01_mapping/).

## Container Types

Containers in this project were created using two different approaches, depending on the tool requirements.

### Single-Tool Containers (Plink, Glimpse, gawk)

**Plink** and **Glimpse** are also used by other stages, so their derivation is documented once in
[`env/shared/README.md`](../shared/README.md) rather than repeated here.

- **Plink**: used for GWAS data format conversion and basic QC filtering.
- **Glimpse**: used for genotype imputation.

**gawk** (`env/01_mapping/gawk.yml` / `gawk.lock`, `bioconda`→`conda-forge::gawk=5.3.1`, frozen to
`community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156`) is only used within
`01_mapping/2_cohort_genotyping`, for text-processing steps between GLIMPSE/BCFtools calls. Same
recovery method as `env/shared/README.md` describes: its `.yml`/`.lock` pair was pulled from the
built image's own `/tmp/conda.yml`/`/tmp/environment.lock`, not reconstructed from guesswork.

### Multi-Tool Containers

Multi-tool containers were created using the [Seqera Wave CLI](https://github.com/seqeralabs/wave-cli) tool, which freezes a complete Conda environment into a Docker image. This approach was used for containers that require multiple interacting tools.

#### Command

```bash
~/Downloads/wave-1.8.1-macos-arm64 --conda-file env/01_mapping/1_High-LowPass_mapping/dog10k_mapping.yml --freeze --await
~/Downloads/wave-1.8.1-macos-arm64 --conda-file env/01_mapping/2_HighPass_SNPfiltering/snpfiltering.yml --freeze --await
```

Frozen image references (as pinned in `code/01_mapping/*/process_definitions.nf`; Wave names each
image after the spec's `name:` field):

| Spec | Image |
|------|-------|
| `1_High-LowPass_mapping/dog10k_mapping.yml` | `community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89` |
| `2_HighPass_SNPfiltering/snpfiltering.yml` | `community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518` |

Both are multi-package builds, so (as `env/shared/README.md` explains) Wave did not bake a
recoverable `environment.lock` into them; re-running the commands above today can resolve newer
transitive dependencies, and the pinned image references are what fix the environment.

#### Parameters

| Parameter | Description |
|-----------|-------------|
| `--conda-file` | Path to the Conda environment YAML file to freeze |
| `--freeze` | Creates a frozen, reproducible container from the environment |
| `--await` | Waits for the container build to complete before returning |

## Environment Files

### [`1_High-LowPass_mapping/dog10k_mapping.yml`](./1_High-LowPass_mapping/dog10k_mapping.yml)

Multi-tool container for the main mapping pipeline. Includes:

| Tool | Version | Purpose |
|------|---------|---------|
| bwa-mem2 | 2.1 | Read alignment |
| gatk4 | 4.4.0.0 | Variant calling |
| picard | 3.1.1 | BAM processing |
| samtools | 1.20 | BAM/BCF manipulation |
| htslib | 1.20 | HTS file I/O |
| python | 3.10.x | Scripting |
| r-base | 4.3.1 | Plotting/analysis |

### [`2_HighPass_SNPfiltering/snpfiltering.yml`](./2_HighPass_SNPfiltering/snpfiltering.yml)

Multi-tool container for high-pass SNP filtering. Includes:

| Tool | Version | Purpose |
|------|---------|---------|
| gatk4 | >=4.3.0.0 | Variant filtering |
| htslib | >=1.19,<1.21 | HTS file I/O |
| bcftools | 1.20 | VCF manipulation |

## Notes

- All containers are built for `linux-64` platform.
- Channels: `conda-forge` and `bioconda`.
- The Seqera Wave CLI (`wave-1.8.1-macos-arm64`) was used on macOS ARM64 (Apple Silicon) to build containers targeting Linux AMD64.
