// ============================================================================
// process_definitions.nf — 3. Merging & filtering
//
// Source: 02_data_processing/3_Merging_filtering/01_gencove_axiom_merging_filtering.sh
//         (Stage A: cohort QC + merge)
// ============================================================================

process QC_FILTER_GENCOVE {
    // QC-filter the LowPass GENCOVE cohort (01_mapping's merged PLINK dataset)
    tag   "gencove"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --keep-allele-order --geno 0.05 --mind 0.05 --maf 0.001 --hwe 'midp' 0.00000000000000000001 --bfile ${prefix} --make-bed --out DA_IMP_GENCOVE_QC
    """

    stub:
    """
    touch DA_IMP_GENCOVE_QC.bed
    touch DA_IMP_GENCOVE_QC.bim
    touch DA_IMP_GENCOVE_QC.fam
    """

    output:
    tuple path("DA_IMP_GENCOVE_QC.bed"), path("DA_IMP_GENCOVE_QC.bim"), path("DA_IMP_GENCOVE_QC.fam"), emit: plink
    tuple val("${task.process}"), val('plink'), eval("plink --version | sed 's/PLINK v//;s/ .*//'"), emit: versions, topic: versions
}


process QC_FILTER_AXIOM {
    // QC-filter the Axiom cohort (2_Axiom_imputation's merged PLINK dataset)
    tag   "axiom"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --allow-no-sex --dog --geno 0.05 --keep-allele-order --mind 0.05 --maf 0.001 --bfile ${prefix} --make-bed --out DA_AFFY_IMP --hwe 'midp' 0.00000000000000000001
    """

    stub:
    """
    touch DA_AFFY_IMP.bed
    touch DA_AFFY_IMP.bim
    touch DA_AFFY_IMP.fam
    """

    output:
    tuple path("DA_AFFY_IMP.bed"), path("DA_AFFY_IMP.bim"), path("DA_AFFY_IMP.fam"), emit: plink
}


process MERGE_DIAGNOSTIC {
    // Diagnostic-only mismatch report before the real merge. `--merge-mode 6` ("no merge, report
    // all mismatching calls") writes only a .diff file — it does not produce .bed/.bim/.fam even
    // with --make-bed specified. The original writes this to the same --out prefix
    // ("DA_MERGED_GENCOVE_AXIOM") that the real merge (--merge-mode 1) below also writes to,
    // meaning the real merge's output overwrites this diagnostic's would-be output — but since
    // mode 6 never actually produces .bed/.bim/.fam, nothing is lost either way. Renamed to its
    // own prefix here purely to avoid a Nextflow output collision with MERGE_COHORTS.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(gencove_bed), path(gencove_bim), path(gencove_fam)
    tuple path(axiom_bed), path(axiom_bim), path(axiom_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def gencove_prefix = gencove_bed.baseName
    def axiom_prefix    = axiom_bed.baseName
    """
    plink --dog --keep-allele-order --bfile ${gencove_prefix} --bmerge ${axiom_prefix} --make-bed --merge-mode 6 --merge-equal-pos --out DA_MERGED_GENCOVE_AXIOM_diagnostic
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM_diagnostic.diff
    """

    output:
    path("DA_MERGED_GENCOVE_AXIOM_diagnostic.diff"), emit: diff
}


process MERGE_COHORTS {
    // The real merge. merge-mode 1 (default): ignore missing calls, set mismatches to missing.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(gencove_bed), path(gencove_bim), path(gencove_fam)
    tuple path(axiom_bed), path(axiom_bim), path(axiom_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def gencove_prefix = gencove_bed.baseName
    def axiom_prefix    = axiom_bed.baseName
    """
    plink --dog --keep-allele-order --bfile ${gencove_prefix} --bmerge ${axiom_prefix} --make-bed --merge-mode 1 --merge-equal-pos --out DA_MERGED_GENCOVE_AXIOM
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM.bed
    touch DA_MERGED_GENCOVE_AXIOM.bim
    touch DA_MERGED_GENCOVE_AXIOM.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM.bed"), path("DA_MERGED_GENCOVE_AXIOM.bim"), path("DA_MERGED_GENCOVE_AXIOM.fam"), emit: plink
}


process QC_MERGED {
    // QC filter on the merged cohort
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --bfile ${prefix} --make-bed --geno 0.05 --maf 0.001 --hwe 'midp' 0.00000000000000000001 --out DA_MERGED_GENCOVE_AXIOM_QC
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM_QC.bed
    touch DA_MERGED_GENCOVE_AXIOM_QC.bim
    touch DA_MERGED_GENCOVE_AXIOM_QC.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM_QC.bed"), path("DA_MERGED_GENCOVE_AXIOM_QC.bim"), path("DA_MERGED_GENCOVE_AXIOM_QC.fam"), emit: plink
}


process FLIP_SCAN_CONVERGENCE {
    // Iterative strand-flip detection and exclusion, converging when a round's flip-scan finds
    // nothing left to exclude — a purely file-based criterion (an empty exclude list) with no
    // human judgment involved, unlike the KING-based duplicate/relatedness exclusion later in
    // this pipeline.
    //
    // Single process with an embedded bash loop, not one Nextflow process per round: DSL2
    // processes can't be invoked in a loop, and there's no fan-out here to lose by keeping it as
    // one task.
    //
    // Two things resolved, not reproduced literally:
    // - Round 1's `--make-pheno DA_IMP_GENCOVE.fam '*'` references a file never created anywhere
    //   in these scripts (only `DA_IMP_GENCOVE_QC.fam` exists, which rounds 2+ correctly use) —
    //   standardized to `DA_IMP_GENCOVE_QC.fam` for every round.
    // - Round 1's `--freq` diagnostic checks on the flipped markers (comparing MAF between the two
    //   platforms) are not reproduced — round 2 has no equivalent, and this is diagnostic-only,
    //   feeding nothing downstream.
    //
    // `.flipscan` has a header row, and `awk '$10 != "NA"'` alone would keep that header line
    // forever (its 10th field is a column label, never literally "NA"), so the derived exclude
    // list would never be truly empty and this loop would never converge. `NR>1` skips the
    // header so the loop converges correctly.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(merged_bed), path(merged_bim), path(merged_fam)
    tuple path(gencove_bed), path(gencove_bim), path(gencove_fam)
    tuple path(axiom_bed), path(axiom_bim), path(axiom_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def merged_prefix  = merged_bed.baseName
    def gencove_prefix = gencove_bed.baseName
    def axiom_prefix   = axiom_bed.baseName
    """
    cur_merged=${merged_prefix}
    cur_gencove=${gencove_prefix}
    cur_axiom=${axiom_prefix}
    round=1

    while :; do
        plink --dog --bfile "\${cur_merged}" --allow-no-sex --make-pheno ${gencove_prefix}.fam '*' --flip-scan --out flipped\${round}
        awk 'NR>1 && \$10 != "NA"' flipped\${round}.flipscan > flipped\${round}.flipscan_noNA
        awk '{print \$2}' flipped\${round}.flipscan_noNA > exclude\${round}_SNPs.txt

        if [ ! -s exclude\${round}_SNPs.txt ]; then
            break
        fi

        plink --dog --bfile "\${cur_gencove}" --exclude exclude\${round}_SNPs.txt --make-bed --out gencove_tmp\${round}
        plink --dog --bfile "\${cur_axiom}" --exclude exclude\${round}_SNPs.txt --make-bed --out axiom_tmp\${round}
        plink --dog --bfile gencove_tmp\${round} --bmerge axiom_tmp\${round} --merge-mode 1 --merge-equal-pos --make-bed --out merged_tmp\${round}

        cur_gencove=gencove_tmp\${round}
        cur_axiom=axiom_tmp\${round}
        cur_merged=merged_tmp\${round}
        round=\$((round + 1))
    done

    cp "\${cur_merged}.bed" DA_MERGED_GENCOVE_AXIOM_QC3.bed
    cp "\${cur_merged}.bim" DA_MERGED_GENCOVE_AXIOM_QC3.bim
    cp "\${cur_merged}.fam" DA_MERGED_GENCOVE_AXIOM_QC3.fam
    """

    stub:
    """
    touch flipped1.flipscan
    touch exclude1_SNPs.txt
    touch DA_MERGED_GENCOVE_AXIOM_QC3.bed
    touch DA_MERGED_GENCOVE_AXIOM_QC3.bim
    touch DA_MERGED_GENCOVE_AXIOM_QC3.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM_QC3.bed"), path("DA_MERGED_GENCOVE_AXIOM_QC3.bim"), path("DA_MERGED_GENCOVE_AXIOM_QC3.fam"), emit: plink
    path("flipped*.flipscan"), emit: flipscan_logs
    path("exclude*_SNPs.txt"), emit: exclude_lists
}


// ============================================================================
// Stage C: sample-level QC — duplicates, depth, KING relatedness, dogID relabeling
// ============================================================================

process QC_FINAL_MERGE {
    // QC filter on the flip-scan-converged merge. Note the --hwe argument order here
    // (p-value then "midp") matches plink's documented `--hwe <p-value> [midp]` form, unlike
    // every other --hwe call in this pipeline (which write "midp" first) — reproduced exactly
    // as each was actually written.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --maf 0.001 --allow-no-sex --geno 0.05 --hwe 0.00000000000000000001 "midp" --bfile ${prefix} --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC4
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM_QC4.bed
    touch DA_MERGED_GENCOVE_AXIOM_QC4.bim
    touch DA_MERGED_GENCOVE_AXIOM_QC4.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM_QC4.bed"), path("DA_MERGED_GENCOVE_AXIOM_QC4.bim"), path("DA_MERGED_GENCOVE_AXIOM_QC4.fam"), emit: plink
}


process BUILD_EXCLUDE_LIST_ROUND1 {
    // Mechanically derive the first-round exclusion list (low depth + known duplicates) from
    // frozen inputs. See bin/build_exclude_list_round1.R.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(fam_qc3modi, arity: '1')
    path(depth_file, arity: '1')
    path(duplicates_file, arity: '1')
    tuple path(qc4_bed), path(qc4_bim), path(qc4_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_exclude_list_round1.R ${fam_qc3modi} ${depth_file} ${duplicates_file} ${qc4_fam}
    """

    stub:
    """
    touch build_exclude_list_round1.log
    touch list_exclude.txt
    """

    output:
    path("list_exclude.txt"),                emit: exclude_list
    path("build_exclude_list_round1.log"),    emit: log
}


process APPLY_EXCLUDE_ROUND1 {
    // Remove low-depth/duplicate samples identified above.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(exclude_list, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --allow-no-sex --remove ${exclude_list} --bfile ${prefix} --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC4B
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM_QC4B.bed
    touch DA_MERGED_GENCOVE_AXIOM_QC4B.bim
    touch DA_MERGED_GENCOVE_AXIOM_QC4B.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM_QC4B.bed"), path("DA_MERGED_GENCOVE_AXIOM_QC4B.bim"), path("DA_MERGED_GENCOVE_AXIOM_QC4B.fam"), emit: plink
}


process PRUNE_FOR_KING {
    // LD-pruned SNP set used by both KING rounds (round 2 reuses this same pruned.prune.in).
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --bfile ${prefix} --indep-pairwise 50 5 0.2 --out pruned
    """

    stub:
    """
    touch pruned.prune.in
    touch pruned.prune.out
    """

    output:
    path("pruned.prune.in"),  emit: prune_in
    path("pruned.prune.out"), emit: prune_out
}


process PCA_ROUND1 {
    // KING-specific PCA/pruned dataset — distinct from the later ALLFAM PCA that feeds
    // 4_Plotting.
    tag   "round1"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(prune_in, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --pca --bfile ${prefix} --allow-no-sex --extract ${prune_in} --out prunedDA_MERGED_GENCOVE_AXIOM_QC4B --make-bed
    """

    stub:
    """
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4B.bed
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4B.bim
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4B.fam
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4B.eigenvec
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4B.eigenval
    """

    output:
    tuple path("prunedDA_MERGED_GENCOVE_AXIOM_QC4B.bed"), path("prunedDA_MERGED_GENCOVE_AXIOM_QC4B.bim"), path("prunedDA_MERGED_GENCOVE_AXIOM_QC4B.fam"), emit: plink
    path("prunedDA_MERGED_GENCOVE_AXIOM_QC4B.eigenvec"), emit: eigenvec
    path("prunedDA_MERGED_GENCOVE_AXIOM_QC4B.eigenval"), emit: eigenval
}


process RUN_KING_ROUND1 {
    // First-round KING relatedness scan — its printed relationship summary is what the original
    // author eyeballed (alongside depth/dogID notes) to manually curate
    // list_exclude_secondexclusion.txt. That curation itself isn't reproducible; this diagnostic
    // run is.
    tag   "round1"
    label 'process_medium'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    king -b ${prefix}.bed --related --rplot --prefix ${prefix} --sexchr 39 > ${prefix}_king.log 2>&1
    """

    stub:
    """
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4B_king.log
    """

    output:
    path("${bed.baseName}_king.log"), emit: king_log
    path("${bed.baseName}*"),         emit: king_files
}


process APPLY_EXCLUDE_ROUND2 {
    // Remove the second round of duplicates — a genuinely manual decision (KING output eyeballed
    // against dog-ID/depth notes, not derivable from any shown code), consumed here as a frozen,
    // versioned input.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(exclude_list_round2, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --allow-no-sex --remove ${exclude_list_round2} --bfile ${prefix} --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC4C
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM_QC4C.bed
    touch DA_MERGED_GENCOVE_AXIOM_QC4C.bim
    touch DA_MERGED_GENCOVE_AXIOM_QC4C.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM_QC4C.bed"), path("DA_MERGED_GENCOVE_AXIOM_QC4C.bim"), path("DA_MERGED_GENCOVE_AXIOM_QC4C.fam"), emit: plink
}


process PCA_ROUND2 {
    // "Double check with KING" — reuses round 1's pruned SNP set, not a fresh pruning.
    tag   "round2"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(prune_in, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --pca --bfile ${prefix} --allow-no-sex --extract ${prune_in} --out prunedDA_MERGED_GENCOVE_AXIOM_QC4C --make-bed
    """

    stub:
    """
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4C.bed
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4C.bim
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4C.fam
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4C.eigenvec
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4C.eigenval
    """

    output:
    tuple path("prunedDA_MERGED_GENCOVE_AXIOM_QC4C.bed"), path("prunedDA_MERGED_GENCOVE_AXIOM_QC4C.bim"), path("prunedDA_MERGED_GENCOVE_AXIOM_QC4C.fam"), emit: plink
    path("prunedDA_MERGED_GENCOVE_AXIOM_QC4C.eigenvec"), emit: eigenvec
    path("prunedDA_MERGED_GENCOVE_AXIOM_QC4C.eigenval"), emit: eigenval
}


process RUN_KING_ROUND2 {
    // Verification re-run confirming round 2's exclusion resolved the relatedness the author was
    // chasing.
    tag   "round2"
    label 'process_medium'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    king -b ${prefix}.bed --related --rplot --prefix ${prefix} --sexchr 39 > ${prefix}_king.log 2>&1
    """

    stub:
    """
    touch prunedDA_MERGED_GENCOVE_AXIOM_QC4C_king.log
    """

    output:
    path("${bed.baseName}_king.log"), emit: king_log
    path("${bed.baseName}*"),         emit: king_files
}


process APPLY_EXCLUDE_LOWPASS12 {
    // Not in the original script. These 12 samples have no Darwin's Ark dogID and no phenotypes;
    // they are not named in any of the three exclude-list files used elsewhere in this pipeline
    // and otherwise ride through QC4C untouched (present in
    // modi_DA_MERGED_GENCOVE_AXIOM_QC4C_forDogIDchange.fam and in the frozen data6 checkpoint as
    // all-NA rows).
    //
    // Removed here from QC4C's own native fam (still keyed by raw sampleID — the dogID relabel
    // hasn't happened yet) via plink's ordinary ID-based --remove, the same mechanism as
    // APPLY_EXCLUDE_ROUND1/2 — not by trimming the dogID-relabel crosswalk fam directly, which
    // would desync its row order from QC4C's .bed/.bim before RELABEL_AND_EXCLUDE_ROUND3's
    // positional --fam substitution (plink matches --fam rows to .bed/.bim rows by position, not
    // by ID).
    //
    // exclude_list_lowpass12 (exclude_12_LowPass_samples.txt) has a header row and a single ID
    // column; reformatted inline to plink's headerless two-column FID/IID --remove format
    // (FID==IID for these samples, same convention list_exclude.txt already uses).
    //
    // A further 7 samples with a similar profile (2 Axiom "-a" secondary-run IDs, 5 GENCOVE SRR
    // accessions) are not named in any provided exclude list and are not removed here.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(exclude_list_lowpass12, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    tail -n +2 ${exclude_list_lowpass12} | awk '{print \$1, \$1}' > exclude_lowpass12_famfmt.txt
    plink --dog --allow-no-sex --remove exclude_lowpass12_famfmt.txt --bfile ${prefix} --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC4C_filtered
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM_QC4C_filtered.bed
    touch DA_MERGED_GENCOVE_AXIOM_QC4C_filtered.bim
    touch DA_MERGED_GENCOVE_AXIOM_QC4C_filtered.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM_QC4C_filtered.bed"), path("DA_MERGED_GENCOVE_AXIOM_QC4C_filtered.bim"), path("DA_MERGED_GENCOVE_AXIOM_QC4C_filtered.fam"), emit: plink
}


process FILTER_DOGID_RELABEL_FAM {
    // Companion to APPLY_EXCLUDE_LOWPASS12: drops the same 12 samples from the dogID-relabel
    // crosswalk fam by ID (column 2, IID), preserving the crosswalk's original row order for
    // everything kept — the same relative order plink's own --remove leaves the genotype side in
    // — so it stays row-for-row aligned with APPLY_EXCLUDE_LOWPASS12's now-filtered QC4C
    // .bed/.bim ahead of RELABEL_AND_EXCLUDE_ROUND3's positional --fam substitution.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    path(dogid_relabel_fam, arity: '1')
    path(exclude_list_lowpass12, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    tail -n +2 ${exclude_list_lowpass12} > exclude_lowpass12_ids.txt
    awk 'NR==FNR{exclude[\$1]; next} !(\$2 in exclude)' exclude_lowpass12_ids.txt ${dogid_relabel_fam} > modi_DA_MERGED_GENCOVE_AXIOM_QC4C_forDogIDchange_filtered.fam
    """

    stub:
    """
    touch modi_DA_MERGED_GENCOVE_AXIOM_QC4C_forDogIDchange_filtered.fam
    """

    output:
    path("modi_DA_MERGED_GENCOVE_AXIOM_QC4C_forDogIDchange_filtered.fam"), emit: fam
}


process RELABEL_AND_EXCLUDE_ROUND3 {
    // Relabel sampleID -> dogID using the frozen, manually-prepared crosswalk fam file, and
    // remove one further duplicate (dog 3094 — "not found among duplicates. Unsure what dog it
    // is. remove:" per the original's own note; a purely manual call, frozen here too).
    //
    // Inputs are the LowPass-12-filtered QC4C and crosswalk fam (see APPLY_EXCLUDE_LOWPASS12/
    // FILTER_DOGID_RELABEL_FAM above) rather than QC4C/the crosswalk fam directly, so the
    // resulting QC5 has 3316 dogs (3328 - 12; none of the 12 overlap the dog-3094 duplicate pair
    // this process also removes), not the original's 3328.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(dogid_relabel_fam, arity: '1')
    path(exclude_list_round3, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --bed ${prefix}.bed --bim ${prefix}.bim --fam ${dogid_relabel_fam} --remove ${exclude_list_round3} --make-bed --out DA_MERGED_GENCOVE_AXIOM_QC5
    """

    stub:
    """
    touch DA_MERGED_GENCOVE_AXIOM_QC5.bed
    touch DA_MERGED_GENCOVE_AXIOM_QC5.bim
    touch DA_MERGED_GENCOVE_AXIOM_QC5.fam
    """

    output:
    tuple path("DA_MERGED_GENCOVE_AXIOM_QC5.bed"), path("DA_MERGED_GENCOVE_AXIOM_QC5.bim"), path("DA_MERGED_GENCOVE_AXIOM_QC5.fam"), emit: plink
}


process BUILD_DATA1 {
    // Join the final QC5 fam with covariates. See bin/build_data1.R.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple path(bed), path(bim), path(fam)
    path(covariates, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_data1.R ${fam} ${covariates}
    """

    stub:
    """
    touch data1.txt
    """

    output:
    path("data1.txt"), emit: data1
}


// ============================================================================
// Stage D: ALLFAM stats, PCA, GRM
// Source: 02_data_processing/3_Merging_filtering/02_create_ALLFAM_check_stats.sh
// ============================================================================

process BUILD_ALLFAM_QC6 {
    // Restrict to the frozen, curated final dog list (all dogs used in any of the GWAS analyses).
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(fam_2584dogs, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --make-bed --maf 0.01 --bfile ${prefix} --keep ${fam_2584dogs} --out ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6
    """

    stub:
    """
    touch ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bed
    touch ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bim
    touch ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam
    """

    output:
    tuple path("ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bed"), path("ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bim"), path("ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam"), emit: plink
}


process PRUNE_ALLFAM {
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --bfile ${prefix} --indep-pairwise 50 5 0.2 --out ALLFAM_pruned
    """

    stub:
    """
    touch ALLFAM_pruned.prune.in
    touch ALLFAM_pruned.prune.out
    """

    output:
    path("ALLFAM_pruned.prune.in"),  emit: prune_in
    path("ALLFAM_pruned.prune.out"), emit: prune_out
}


process PCA_ALLFAM {
    // This is the PCA that actually feeds 4_Plotting/plot_PCA.R — distinct from Stage C's
    // KING-specific PCA rounds.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple path(bed), path(bim), path(fam)
    path(prune_in, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink --dog --pca --bfile ${prefix} --allow-no-sex --extract ${prune_in} --out prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 --make-bed
    """

    stub:
    """
    touch prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bed
    touch prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bim
    touch prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam
    touch prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.eigenvec
    touch prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.eigenval
    """

    output:
    tuple path("prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bed"), path("prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.bim"), path("prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.fam"), emit: plink
    path("prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.eigenvec"), emit: eigenvec
    path("prunedALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.eigenval"), emit: eigenval
}


process BUILD_GRM {
    tag   "cohort"
    label 'process_high'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple path(bed), path(bim), path(fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    gcta64 --bfile ${prefix} --autosome --make-grm --out ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6 --autosome-num 38
    """

    stub:
    """
    touch ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.grm.bin
    touch ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.grm.N.bin
    touch ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.grm.id
    """

    output:
    path("ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.grm.bin"),   emit: grm_bin
    path("ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.grm.N.bin"), emit: grm_n_bin
    path("ALLFAM_DA_MERGED_GENCOVE_AXIOM_QC6.grm.id"),    emit: grm_id
}


process DEPTH_HISTOGRAMS {
    // Diagnostic depth-distribution plots, dated later than the rest of this script in the
    // original. Feeds nothing downstream. See bin/depth_histograms.R.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(fam_qc3modi, arity: '1')
    path(depth_file, arity: '1')
    tuple path(allfam_bed), path(allfam_bim), path(allfam_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    depth_histograms.R ${fam_qc3modi} ${depth_file} ${allfam_fam}
    """

    stub:
    """
    touch depth_histograms.log
    touch depth_histograms.png
    """

    output:
    path("depth_histograms.log"), emit: log
    path("depth_histograms.png"), emit: plot
}


// ============================================================================
// SIZE — control-phenotype fam file, run via mlma alongside F1/F2/F3. Not part of the original
// 01/02/03/04 shell/R sources for this pipeline — see data/04_gwas/getSIZE_pheno_genofiles.r.
// That file's own second section (STUCK, a different Darwin's Ark item) is not reproduced here.
// ============================================================================

process BUILD_SIZE_FAM {
    // See bin/build_size_fam.R for the write-then-read-back idiom this resolves (reuses this
    // pipeline's own data1.txt rather than round-tripping through disk).
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(data1, arity: '1')
    path(size_pheno_file, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_size_fam.R ${data1} ${size_pheno_file}
    """

    stub:
    """
    touch SIZE_DA_MERGED_GENCOVE_AXIOM_QC5.fam
    touch build_size_fam.log
    """

    output:
    path("SIZE_DA_MERGED_GENCOVE_AXIOM_QC5.fam"), emit: fam
    path("build_size_fam.log"), emit: log
}


process FILTER_SIZE_QC6 {
    // Forked from FILTER_FACTOR_QC6 rather than reused directly: that process's output filenames
    // hardcode a "CCD" prefix (CCD${factor_id}_...), but the real, already-committed
    // data/04_gwas/SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam has no such prefix — matched here exactly.
    // No source script documents SIZE's own QC6 plink invocation (getSIZE_pheno_genofiles.r only
    // builds the fam); this reuses FILTER_FACTOR_QC6's exact command, since README_04.gwas.md
    // groups SIZE under the same "1_mlma-loco (factors)" GWAS model as F1/F2/F3.
    tag   "SIZE"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    path(size_fam, arity: '1')
    tuple path(qc5_bed), path(qc5_bim), path(qc5_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plink --dog --make-bed --prune --maf 0.01 --bed ${qc5_bed} --bim ${qc5_bim} --fam ${size_fam} --out SIZE_DA_MERGED_GENCOVE_AXIOM_QC6
    awk '{print \$1,\$2,\$6}' SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam > SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.phen
    """

    stub:
    """
    touch SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.bed
    touch SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.bim
    touch SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam
    touch SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.phen
    """

    output:
    tuple path("SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.bed"), path("SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.bim"), path("SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.fam"), emit: plink
    path("SIZE_DA_MERGED_GENCOVE_AXIOM_QC6.phen"), emit: phen
}


process BUILD_SIZE_COVARIATES {
    // Forked from BUILD_FACTOR_COVARIATES rather than reused directly: that process's output
    // filenames hardcode a "CCD" prefix, but SIZE's real files don't use one (same reasoning as
    // FILTER_SIZE_QC6). Reuses bin/build_factor_covariates.R unchanged — that script takes its
    // output prefix as a plain argument, already fully generic. SIZE's age/sex covariates are
    // built the same way as the factors' covariates.
    tag   "SIZE"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple path(qc6_bed), path(qc6_bim), path(qc6_fam)
    path(data1, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_factor_covariates.R ${qc6_fam} ${data1} SIZE_DA_MERGED_GENCOVE_AXIOM_QC6
    """

    stub:
    """
    touch ageSIZE_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar
    touch sexSIZE_DA_MERGED_GENCOVE_AXIOM_QC6.covar
    touch build_SIZE_DA_MERGED_GENCOVE_AXIOM_QC6_covariates.log
    """

    output:
    path("ageSIZE_DA_MERGED_GENCOVE_AXIOM_QC6.qcovar"), emit: qcovar
    path("sexSIZE_DA_MERGED_GENCOVE_AXIOM_QC6.covar"),  emit: covar
    path("build_SIZE_DA_MERGED_GENCOVE_AXIOM_QC6_covariates.log"), emit: log
}


// ============================================================================
// Stage E: per-factor phenotype files (F1/F2/F3)
// Source: 02_data_processing/3_Merging_filtering/03_create_files_per_factor_mlma.R
//
// `data6_2024-10-14.txt` is consumed here as a frozen input rather than recomputed from data1 +
// F1/F2/F3_CCD3F.txt (the join the original script itself performs immediately before reading
// this file back in, overwriting the freshly-computed result). Recomputing is not done here: the
// original's `age = data5$age.x, sex = data5$sex_numeric.x` column references imply a name
// collision between data1 and the per-factor CCD3F tables that doesn't reproduce with the CCD3F
// files this pipeline's own 1_Survey_data stage produces (which carry no `age`/`sex_numeric`
// columns to collide with) — recomputing would silently drop those columns or error on first
// downstream use, and there is no committed historical F1/F2/F3 CCD3F file to check the original
// collision against.
//
// data6 also needs restricting to QC5's own dogID set/order before BUILD_FACTOR_FAM/
// BUILD_ITEM_FAM can use it — see FILTER_DATA6_TO_QC5's header for why (a consequence of the
// LowPass-12 exclusion in Stage C). BUILD_ITEM_GRAB (Stage F) keeps using the *unfiltered* data6,
// since its join there is ID-based, not a positional --fam substitution.
// ============================================================================

process FILTER_DATA6_TO_QC5 {
    // Restricts + reorders the frozen data6 checkpoint to exactly QC5's own dogID set, in QC5's
    // own row order. QC5 is 3316 rows after Stage C's LowPass-12 exclusion, while data6 is still
    // 3328 (frozen, untouched by that exclusion), and BUILD_FACTOR_FAM/BUILD_ITEM_FAM's output
    // feeds a strictly positional plink --fam substitution downstream (FILTER_FACTOR_QC6/
    // FILTER_ITEM_QC6) that requires exact row-count and row-order agreement with QC5's
    // .bed/.bim.
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    path(data6, arity: '1')
    tuple path(qc5_bed), path(qc5_bim), path(qc5_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    filter_data6_to_qc5.R ${data6} ${qc5_fam}
    """

    stub:
    """
    touch data6_filtered.txt
    """

    output:
    path("data6_filtered.txt"), emit: data6
}


process BUILD_FACTOR_FAM {
    // Per-factor NA-count-filtered phenotype fam file. See bin/build_factor_fam.R for the
    // na_count bug fix (bare `na_count` -> `na_count_F1`) applied for F1.
    tag   "${factor_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple val(factor_id), val(item_cols), val(na_threshold)
    path(data6, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_factor_fam.R ${data6} ${factor_id} '${item_cols}' ${na_threshold}
    """

    stub:
    """
    touch CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC5.fam
    touch build_${factor_id}_fam.log
    touch ${factor_id}_na_counts.png
    """

    output:
    tuple val(factor_id), path("CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC5.fam"), emit: fam
    path("build_${factor_id}_fam.log"), emit: log
    path("${factor_id}_na_counts.png"), emit: na_counts_plot
}


process FILTER_FACTOR_QC6 {
    // Note: unlike the per-item filter below, this one has no --allow-no-sex — reproduced
    // exactly as each was actually written.
    tag   "${factor_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(factor_id), path(factor_fam)
    tuple path(qc5_bed), path(qc5_bim), path(qc5_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plink --dog --make-bed --prune --maf 0.01 --bed ${qc5_bed} --bim ${qc5_bim} --fam ${factor_fam} --out CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6
    awk '{print \$1,\$2,\$6}' CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.fam > CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.phen
    """

    stub:
    """
    touch CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bed
    touch CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bim
    touch CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.fam
    touch CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.phen
    """

    output:
    tuple val(factor_id), path("CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bed"), path("CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bim"), path("CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.fam"), emit: plink
    path("CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6.phen"), emit: phen
}


process BUILD_FACTOR_COVARIATES {
    // Joins the QC6 fam with data1 (see bin/build_factor_covariates.R for the data4 -> data1
    // fix) to build age qcovar / sex covar files feeding 04_gwas (not yet converted).
    tag   "${factor_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple val(factor_id), path(qc6_bed), path(qc6_bim), path(qc6_fam)
    path(data1, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = "CCD${factor_id}_DA_MERGED_GENCOVE_AXIOM_QC6"
    """
    build_factor_covariates.R ${qc6_fam} ${data1} ${prefix}
    """

    stub:
    """
    touch age${factor_id}_placeholder.qcovar
    touch sex${factor_id}_placeholder.covar
    touch build_placeholder_covariates.log
    """

    output:
    path("age*.qcovar"), emit: qcovar
    path("sex*.covar"),  emit: covar
    path("build_*_covariates.log"), emit: log
}


// ============================================================================
// Stage F: per-item phenotype files (14 items)
// Source: 02_data_processing/3_Merging_filtering/04_create_files_per_item_polmm.sh
//
// The original shows one worked example (item155) with a "repeat for each item" comment —
// genericized here across all 14 items via the item_id parameter.
// ============================================================================

process BUILD_ITEM_FAM {
    tag   "item${item_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    val(item_id)
    path(data6, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_item_fam.R ${data6} ${item_id}
    """

    stub:
    """
    touch CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC5.fam
    touch build_item${item_id}_fam.log
    """

    output:
    tuple val(item_id), path("CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC5.fam"), emit: fam
    path("build_item${item_id}_fam.log"), emit: log
}


process FILTER_ITEM_QC6 {
    tag   "item${item_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(item_id), path(item_fam)
    tuple path(qc5_bed), path(qc5_bim), path(qc5_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plink --dog --allow-no-sex --make-bed --prune --maf 0.01 --bed ${qc5_bed} --bim ${qc5_bim} --fam ${item_fam} --out CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6
    """

    stub:
    """
    touch CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bed
    touch CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bim
    touch CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.fam
    """

    output:
    tuple val(item_id), path("CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bed"), path("CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.bim"), path("CCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.fam"), emit: plink
}


process PRUNE_ITEM_LD {
    tag   "item${item_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(item_id), path(qc6_bed), path(qc6_bim), path(qc6_fam)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = qc6_bed.baseName
    """
    plink --bfile ${prefix} --dog --indep-pairwise 50 5 0.2 --out CCDitem${item_id}_indepSNP
    """

    stub:
    """
    touch CCDitem${item_id}_indepSNP.prune.in
    touch CCDitem${item_id}_indepSNP.prune.out
    """

    output:
    tuple val(item_id), path("CCDitem${item_id}_indepSNP.prune.in"), emit: prune_in
    path("CCDitem${item_id}_indepSNP.prune.out"), emit: prune_out
}


process BUILD_ITEM_EIGENINPUT {
    tag   "item${item_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(item_id), path(qc6_bed), path(qc6_bim), path(qc6_fam), path(prune_in)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = qc6_bed.baseName
    """
    plink --bfile ${prefix} --dog --extract ${prune_in} --make-bed --out CCDitem${item_id}_EigenInput
    """

    stub:
    """
    touch CCDitem${item_id}_EigenInput.bed
    touch CCDitem${item_id}_EigenInput.bim
    touch CCDitem${item_id}_EigenInput.fam
    """

    output:
    tuple val(item_id), path("CCDitem${item_id}_EigenInput.bed"), path("CCDitem${item_id}_EigenInput.bim"), path("CCDitem${item_id}_EigenInput.fam"), emit: plink
}


process BUILD_ITEM_GRAB {
    // Combined phenotype/covariate file feeding POLMM. See bin/build_item_grab.R for the
    // sex.x/sex.y join-collision note.
    tag   "item${item_id}"
    label 'process_single'
    container 'community.wave.seqera.io/library/merging_filtering_tools:0e56b6d8b7fd7863'

    input:
    tuple val(item_id), path(qc6_bed), path(qc6_bim), path(qc6_fam)
    path(data6, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    build_item_grab.R ${qc6_fam} ${data6} ${item_id}
    """

    stub:
    """
    touch grabCCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.txt
    touch build_item${item_id}_grab.log
    """

    output:
    path("grabCCDitem${item_id}_DA_MERGED_GENCOVE_AXIOM_QC6.txt"), emit: grab
    path("build_item${item_id}_grab.log"), emit: log
}


// ============================================================================
// Stage G: cross-phenotype dog/SNP-count check — a barrier over all 17 (3 factor + 14 item)
// QC6 datasets.
// Source: 02_data_processing/3_Merging_filtering/02_create_ALLFAM_check_stats.sh, lines 26-84
// ============================================================================

process CHECK_PHENOTYPE_COUNTS {
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    val(phenotype_ids)
    path(fams, arity: '1..*')
    path(bims, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    def ids = phenotype_ids.join(' ')
    """
    check_phenotype_counts.sh ${fams} -- ${bims} -- ${ids} > phenotype_counts.log
    """

    stub:
    """
    touch phenotype_counts.log
    """

    output:
    path("phenotype_counts.log"), emit: log
}
