// ============================================================================
// process_definitions.nf — 04. GWAS — 3. Finemap (SuSiE, no functional priors)
//
// Source: 04_gwas/3_finemap_susie/{1_finemap_create_genofiles.sh,
//          2a_run_harmonize_sumstats_CCD.sbatch, harmonize_sumstats_CCD_SIZE.py,
//          2b_munge_sumstats_parquet.sh, finemapNOFUNCT_multiple_CCDregions{1,2,3,4}.sh}
//
// Structurally different from 1_mlma-loco/2_polmm: instead of one worked example genericized
// across a phenotype/item list, the finemapper.py calls are ~100 individual, hand-picked
// (phenotype, chr, start, end) windows accumulated over an iterative, exploratory process — some
// loci re-run at 2-4 progressively narrower widths. Every one of these calls is reproduced
// literally (not just the apparently-final/narrowest one per locus). The full table was extracted
// mechanically (a one-off parsing script, not hand-typed) from the 4 original shell scripts into
// assets/finemap_regions.csv (104 rows) to avoid transcription error — see that file's own header.
//
// STUCK is excluded from harmonize/munge/finemap entirely, matching the literal original: no
// finemap-region windows exist for it anywhere, and it was excluded from GWAS entirely per
// 3_Merging_filtering's own scope. Its deposited fam file (STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.fam,
// referenced below) is used only as a sample-inclusion mask for the one shared genofile all
// phenotypes finemap against — not as an analysed phenotype.
//
// SIZE, by contrast, IS harmonized/munged/finemapped here, despite
// 2a_run_harmonize_sumstats_CCD.sbatch's own SIZE line being commented out in the original. No
// finemap-region windows for SIZE exist in any of the 4 region-table scripts either — a finemap
// window is a human-picked choice per locus, so there's no generic "same shape as the factors"
// command for it the way there was for clumping. data/04_gwas/SIZE_260603_CSnonfunct.xlsx is
// direct evidence SIZE was finemapped in the real analysis — its own already-run region filenames
// double as the missing window list, mechanically extracted into assets/finemap_regions_size.csv
// (9 rows) by extract_size_regions.py. The credible-set boundary columns in that spreadsheet are
// NOT used as the finemap window — they're the credible set's own narrower bounds, not the region
// analysed to find it (see extract_size_regions.py's header).
//
// New container required: PolyFun (github.com/omerwe/polyfun) is an external tool with no
// versioned releases (rolling git history only), not vendored in this repo, not pip/conda-
// installable as a package. Pinned to commit 3e657a1066 (2024-07-28), built from PolyFun's own
// polyfun.yml (conda-forge only, python=3.8, r-susier==0.11.92 via rpy2) plus a `git clone` +
// `git checkout` layered on top (see env/04_gwas/3_finemap_susie/polyfun.yml and
// env/04_gwas/README.md). Reused for harmonize/munge/finemap alike, since harmonize only needs
// pandas (already in PolyFun's own env).
//
// The genofile-creation step reconstructs the combined fam file using 3_Merging_filtering's own
// published QC5 fam (DA_MERGED_GENCOVE_AXIOM_QC5.fam) — see bin/build_combined_fam.R's header.
//
// The original's `--hwe 'midp' 0.00000000000000000001` argument order (modifier before the
// p-value) is reproduced literally below despite plink 1.9's documented syntax being
// `--hwe <p-value> [midp]` (value first); it is unverified whether plink accepts this reversed
// order or silently misparses it.
// ============================================================================

process BUILD_COMBINED_FAM {
    // See bin/build_combined_fam.R's header for the famQC5 reconstruction and the
    // sample-inclusion-mask (not phenotype) role of STUCK here.
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(qc5_fam, arity: '1')
    path(size_fam, arity: '1')
    path(stuck_fam, arity: '1')
    path(allfam_fam, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_combined_fam.R ${qc5_fam} ${size_fam} ${stuck_fam} ${allfam_fam} CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC5.fam
    """

    stub:
    """
    touch CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC5.fam
    """

    output:
    path("CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC5.fam"), emit: fam
}


process FILTER_COMBINED_QC6 {
    label 'process_high'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(qc5_bed), path(qc5_bim)
    path(combined_fam, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plink --dog --make-bed --allow-no-sex --prune --geno 0.05 --mind 0.05 --maf 0.01 --hwe 'midp' 0.00000000000000000001 --bed ${qc5_bed} --bim ${qc5_bim} --fam ${combined_fam} --out CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6
    """

    stub:
    """
    touch CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.bed
    touch CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.bim
    touch CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.fam
    touch CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.log
    """

    output:
    tuple path("CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.bed"), path("CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.bim"), path("CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.fam"), emit: plink
    path("CCD_SIZE_STUCK_DA_MERGED_GENCOVE_AXIOM_QC6.log"), emit: log
}


process SPLIT_GENOFILE_BY_CHR {
    // Fans out over chr 1-38 (dog autosomes) — matches the original's own `for chr in {1..38}` loop.
    tag   "chr${chr}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(chr), path(qc6_bed), path(qc6_bim), path(qc6_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = qc6_bed.baseName
    """
    plink --dog --bfile ${prefix} --chr ${chr} --make-bed --out ${prefix}.${chr}
    """

    stub:
    """
    touch ${qc6_bed.baseName}.${chr}.bed
    touch ${qc6_bed.baseName}.${chr}.bim
    touch ${qc6_bed.baseName}.${chr}.fam
    """

    output:
    tuple val(chr), path("${qc6_bed.baseName}.${chr}.bed"), path("${qc6_bed.baseName}.${chr}.bim"), path("${qc6_bed.baseName}.${chr}.fam"), emit: plink
}


process HARMONIZE_SUMSTATS {
    // All 18 GWAS'd phenotypes (3 factors + 14 items + SIZE) — STUCK excluded, see this file's header.
    tag   "${pheno_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/polyfun:cb5ec2572b0e5a72'

    input:
    tuple val(pheno_id), path(sumstats)
    path(combined_bim, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    harmonize_sumstats.py ${sumstats} ${pheno_id}_sumstats_harmonized.txt ${combined_bim} ${pheno_id}_swapped_snps.txt
    """

    stub:
    """
    touch ${pheno_id}_sumstats_harmonized.txt
    touch ${pheno_id}_swapped_snps.txt
    """

    output:
    tuple val(pheno_id), path("${pheno_id}_sumstats_harmonized.txt"), emit: harmonized
    path("${pheno_id}_swapped_snps.txt"), emit: swapped_snps
}


process MUNGE_SUMSTATS {
    tag   "${pheno_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/polyfun:cb5ec2572b0e5a72'

    input:
    tuple val(pheno_id), path(harmonized_sumstats)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 /opt/polyfun/munge_polyfun_sumstats.py --sumstats ${harmonized_sumstats} --keep-hla --min-info 0 --out ${pheno_id}_sumstats_munged.parquet
    """

    stub:
    """
    touch ${pheno_id}_sumstats_munged.parquet
    """

    output:
    tuple val(pheno_id), path("${pheno_id}_sumstats_munged.parquet"), emit: munged
}


process FINEMAP_REGION {
    // One task per row of assets/finemap_regions.csv (104 rows, mechanically extracted from the
    // original's 4 finemapNOFUNCT_multiple_CCDregions*.sh scripts) + assets/finemap_regions_size.csv
    // (9 SIZE rows, mechanically extracted from SIZE_260603_CSnonfunct.xlsx) — see this file's header.
    tag   "${row.out_base}"
    label 'process_high'
    container 'community.wave.seqera.io/library/polyfun:cb5ec2572b0e5a72'

    input:
    tuple val(row), path(chr_bed), path(chr_bim), path(chr_fam), path(munged_parquet)

    when:
    task.ext.when == null || task.ext.when

    script:
    def geno_prefix = chr_bed.baseName
    """
    python3 /opt/polyfun/finemapper.py --geno ${geno_prefix} --sumstats ${munged_parquet} --non-funct --n ${row.n} --chr ${row.chr} --start ${row.start} --end ${row.end} --method susie --allow-missing --max-num-causal 5 --out ${row.out_base}
    """

    stub:
    """
    touch ${row.out_base}
    """

    output:
    path("${row.out_base}"), emit: result
}
