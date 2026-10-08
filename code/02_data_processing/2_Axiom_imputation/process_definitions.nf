// ============================================================================
// process_definitions.nf — 2. Axiom imputation
//
// Source scripts:
//   02_data_processing/2_Axiom_imputation/01_axiom_liftover.sh    (LIFTOVER_TO_CANFAM4)
//   02_data_processing/2_Axiom_imputation/02_axiom_impute.sh      (SPLIT_REF_PANEL_BY_CHR,
//                                                                   RECODE_AND_RENAME_CHR_VCF,
//                                                                   CONFORM_GT, BEAGLE_IMPUTE,
//                                                                   QC_AND_CONVERT_CHR)
//   02_data_processing/2_Axiom_imputation/03_axiomimp_to_plink.sh (QC_AND_CONVERT_CHR,
//                                                                   MERGE_CHR_PLINK_AXIOM)
//
// Process granularity: each process below corresponds to a group of adjacent original commands
// that always run together (same scope, no fan-out or other consumer between them) and share a
// similar resource profile. Steps kept separate — SPLIT_REF_PANEL_BY_CHR, CONFORM_GT,
// BEAGLE_IMPUTE, MERGE_CHR_PLINK_AXIOM — either have a genuinely different resource profile
// (BEAGLE_IMPUTE needs process_high; the steps around it don't) or are independently reusable/
// retry-prone enough that merging them would cost more in lost -resume granularity than it would
// save in container-start overhead.
// ============================================================================

process LIFTOVER_TO_CANFAM4 {
    // Full canFam3 -> canFam4 liftover chain for the Axiom array genotypes.
    // Source: 02_data_processing/2_Axiom_imputation/01_axiom_liftover.sh
    //
    // Merged from 11 adjacent steps (filter to 411 samples, sex-check diagnostic, build liftover
    // BED, liftOver, derive marker bookkeeping files, filter+flip lifted markers, marker-quality
    // diagnostics, exclude mismatched-chr markers, update map/chr, re-sort positions, rename SNP
    // IDs): a strictly linear, cohort-level, single-plinkset chain with no fan-out and no other
    // consumer of any intermediate file.
    //
    // - `.bim` is used in place of a `.map` file the awk command reads (`affy_round3.map`) but
    //   which `plink --make-bed` never produces; `.bim`'s chr/id/cM/pos columns 1-4 match.
    // - The marker-quality check `$5 == !"A/C/T/G"` is preserved verbatim; awk evaluates this as
    //   `$5 == 0`, duplicating the A1_is_0 check rather than testing for non-ACGT alleles as its
    //   filename implies.
    // - `--hwe 'midp' <p>`: argument order preserved as written (see the pipeline README).
    // - The sex check and marker-quality checks are diagnostics whose output nothing downstream
    //   reads — reproduced anyway.
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e'

    input:
    tuple path(bed), path(bim), path(fam)
    path(ind_to_keep, arity: '1')
    path(chain, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    ## 1. Filter to the 411 kept samples
    plink --bfile ${prefix} --dog --make-bed --keep ${ind_to_keep} --out affy_round3

    ## 2. Sex check (diagnostic, unused downstream)
    plink --dog --bfile affy_round3 --split-x 6600000 123798852 --indep-pairphase 20000 2000 0.5 --check-sex 0.5 0.8 --out sex_affy_round3

    ## 3. Build liftover BED
    awk '{print "chr"\$1,\$4-1,\$4,\$2,"0","+"}' affy_round3.bim > affy_round3_cf3_BED.bed

    ## 4. liftOver: canFam3 -> canFam4
    liftOver affy_round3_cf3_BED.bed ${chain} output.bed unlifted.bed

    ## 5. Derive marker bookkeeping files
    awk '{print \$4}' output.bed > lifted_markers.txt

    grep - output.bed > flip_markers.txt
    awk '{print \$4}' flip_markers.txt > flip_markers_list.txt

    awk -F'\\t' '{split(\$4, a, ":"); if (\$1 != a[1]) print}' output.bed > no_match.txt
    awk '{print \$4}' no_match.txt > no_match_markers_remove.txt

    awk '{print \$4,\$3}' output.bed > map_canfam4_build.txt
    awk '{print \$4,\$1}' output.bed > chr_canfam4_build.txt

    ## 6. Filter to lifted, ACGT-only SNPs; flip strand
    plink --bfile affy_round3 \\
        --extract lifted_markers.txt \\
        --flip flip_markers_list.txt \\
        --geno 0.05 \\
        --maf 0.01 \\
        --hwe 'midp' 0.000000000000001 \\
        --make-bed \\
        --dog \\
        --out affy_cf3 \\
        --snps-only 'just-acgt'

    ## 7. Marker-quality checks (diagnostic, unused downstream)
    cat affy_cf3.bim | awk '{if (\$5 == "0") print \$0;}' > A1_is_0.txt
    cat affy_cf3.bim | awk '{if (\$6 == "0") print \$0;}' > A2_is_0.txt
    cat affy_cf3.bim | awk '{if (\$5 == !"A/C/T/G") print \$0;}' > A1_notACTG.txt
    awk 'length(\$5) > 1 || length(\$6) > 1' affy_cf3.bim > indels.txt

    ## 8. Exclude markers lifted to an unexpected chromosome
    plink --bfile affy_cf3 --exclude no_match_markers_remove.txt --dog --make-bed --out affy_cf3_a

    ## 9. Move positions/chromosomes onto canFam4
    plink --bfile affy_cf3_a --update-map map_canfam4_build.txt --update-chr chr_canfam4_build.txt --dog --make-bed --out affy_cf4

    ## 10. Re-sort variants by position
    plink --bfile affy_cf4 --dog --make-bed --out affy_cf4b

    ## 11. Rename variant IDs to chr:pos (final DA_AFFY plinkset)
    awk '{print \$2,"chr"\$1":"\$4}' affy_cf4b.bim > newSNPID.txt
    plink --bfile affy_cf4b --update-name newSNPID.txt --dog --make-bed --out DA_AFFY
    """

    stub:
    """
    touch sex_affy_round3.sexcheck
    touch A1_is_0.txt
    touch A2_is_0.txt
    touch A1_notACTG.txt
    touch indels.txt
    touch DA_AFFY.bed
    touch DA_AFFY.bim
    touch DA_AFFY.fam
    """

    output:
    tuple path("DA_AFFY.bed"), path("DA_AFFY.bim"), path("DA_AFFY.fam"), emit: plink
    path("sex_affy_round3.sexcheck"),                                   emit: sexcheck
    path("A1_is_0.txt"),                                                emit: a1_is_0
    path("A2_is_0.txt"),                                                emit: a2_is_0
    path("A1_notACTG.txt"),                                             emit: a1_not_actg
    path("indels.txt"),                                                 emit: indels
    tuple val("${task.process}"), val('plink'),    eval("plink --version | sed 's/PLINK v//;s/ .*//'"), emit: versions_plink,    topic: versions
    tuple val("${task.process}"), val('liftOver'), eval("liftOver 2>&1 | sed -n '1p'"),                  emit: versions_liftover, topic: versions
}


process SPLIT_REF_PANEL_BY_CHR {
    // Materialize one chromosome slice of the Dog10K phased reference panel.
    // Source comment: "make one ref per chr: DO THIS ONCE" (02_axiom_impute.sh)
    tag   "${chr}"
    label 'process_low'
    container 'community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e'

    input:
    tuple val(chr), val(start), val(stop)
    tuple path(ref_bcf), path(ref_csi)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view ${ref_bcf} -r ${chr} -Oz -o modi_AutoAndXPAR.Dog10K.Phased_${chr}.vcf.gz
    bcftools index modi_AutoAndXPAR.Dog10K.Phased_${chr}.vcf.gz
    """

    stub:
    """
    touch modi_AutoAndXPAR.Dog10K.Phased_${chr}.vcf.gz
    touch modi_AutoAndXPAR.Dog10K.Phased_${chr}.vcf.gz.csi
    """

    output:
    tuple val(chr), val(start), val(stop),
          path("modi_AutoAndXPAR.Dog10K.Phased_${chr}.vcf.gz"),
          path("modi_AutoAndXPAR.Dog10K.Phased_${chr}.vcf.gz.csi"),                                        emit: ref_chr
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions, topic: versions
}


process RECODE_AND_RENAME_CHR_VCF {
    // Recode one chromosome of DA_AFFY to VCF, orienting alleles to the reference FASTA, then
    // rename its bare integer chromosome code (e.g. "1") to the "chr1" form used by the
    // reference panel.
    // Source: 02_data_processing/2_Axiom_imputation/02_axiom_impute.sh
    //
    // Merged from 2 originally-separate steps: both lightweight, always run back-to-back per
    // chromosome, with no other consumer of the intermediate plain VCF.
    //
    // The original also runs `bcftools view DA_AFFY_$chrN.vcf -Oz -o DA_AFFY_$chrN.vcf.gz` (+
    // index) between these two steps, then never reads that compressed copy again — the rename
    // command re-reads the plain .vcf directly. Reproduced anyway as an unused leaf output.
    tag   "${chr}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e'

    input:
    tuple val(chr), val(start), val(stop)
    tuple path(bed), path(bim), path(fam)
    tuple path(assembly_ref), path(assembly_ref_fai)
    path(chr_rename_map, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = bed.baseName
    """
    plink2 --bfile ${prefix} --dog --chr ${chr} --recode vcf id-paste=iid --fa ${assembly_ref} --ref-from-fa --out DA_AFFY_${chr}

    bcftools view DA_AFFY_${chr}.vcf -Oz -o DA_AFFY_${chr}.vcf.gz
    bcftools index DA_AFFY_${chr}.vcf.gz

    bcftools annotate --rename-chrs ${chr_rename_map} DA_AFFY_${chr}.vcf -Oz -o DA_AFFY_${chr}_rename.vcf.gz
    bcftools index DA_AFFY_${chr}_rename.vcf.gz
    """

    stub:
    """
    touch DA_AFFY_${chr}.vcf.gz
    touch DA_AFFY_${chr}.vcf.gz.csi
    touch DA_AFFY_${chr}_rename.vcf.gz
    touch DA_AFFY_${chr}_rename.vcf.gz.csi
    """

    output:
    tuple val(chr), val(start), val(stop),
          path("DA_AFFY_${chr}_rename.vcf.gz"), path("DA_AFFY_${chr}_rename.vcf.gz.csi"),                 emit: vcf
    tuple path("DA_AFFY_${chr}.vcf.gz"), path("DA_AFFY_${chr}.vcf.gz.csi"),                               emit: unused_plain_vcf
    tuple val("${task.process}"), val('plink2'), eval("plink2 --version | sed 's/.*PLINK v//;s/ .*//'"), emit: versions, topic: versions
}


process CONFORM_GT {
    // Reconcile alleles/strand against the reference panel before phasing/imputation
    tag   "${chr}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e'

    input:
    tuple val(chr), val(start), val(stop), path(gt_vcf), path(gt_csi), path(ref_vcf), path(ref_csi)
    path(conform_gt_jar, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    java -jar ${conform_gt_jar} \\
        ref=${ref_vcf} \\
        gt=${gt_vcf} \\
        chrom=${chr} \\
        match=POS \\
        out=DA_AFFY_conf_${chr} \\
        strict=TRUE
    bcftools index DA_AFFY_conf_${chr}.vcf.gz
    """

    stub:
    """
    touch DA_AFFY_conf_${chr}.vcf.gz
    touch DA_AFFY_conf_${chr}.vcf.gz.csi
    """

    output:
    tuple val(chr), val(start), val(stop), path("DA_AFFY_conf_${chr}.vcf.gz"), path("DA_AFFY_conf_${chr}.vcf.gz.csi"), emit: vcf
    tuple val("${task.process}"), val('openjdk'), eval("java -version 2>&1 | sed -n '1p'"),                            emit: versions, topic: versions
}


process BEAGLE_IMPUTE {
    // Phase and impute against the Dog10K reference panel
    tag   "${chr}"
    label 'process_high'
    container 'community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e'

    input:
    tuple val(chr), val(start), val(stop), path(conf_vcf), path(conf_csi), path(ref_vcf), path(ref_csi), path(genetic_map)
    path(beagle_jar, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    java -Xmx50g -jar ${beagle_jar} \\
        gt=${conf_vcf} \\
        chrom=${chr}:${start}-${stop} \\
        map=${genetic_map} \\
        out=DA_AFFYimp_${chr} \\
        ref=${ref_vcf} \\
        iterations=40 \\
        burnin=30
    bcftools index DA_AFFYimp_${chr}.vcf.gz
    """

    stub:
    """
    touch DA_AFFYimp_${chr}.vcf.gz
    touch DA_AFFYimp_${chr}.vcf.gz.csi
    """

    output:
    tuple val(chr), path("DA_AFFYimp_${chr}.vcf.gz"), path("DA_AFFYimp_${chr}.vcf.gz.csi"), emit: vcf
}


process QC_AND_CONVERT_CHR {
    // QC-filter on imputation quality, rename variant IDs, run diagnostic counts, and convert
    // one chromosome's imputed VCF to PLINK binary format.
    // Source: 02_data_processing/2_Axiom_imputation/02_axiom_impute.sh (QC/rename/counts) and
    //         03_axiomimp_to_plink.sh (VCF -> PLINK)
    //
    // Merged from 4 originally-separate steps: all lightweight, always run back-to-back per
    // chromosome, similar resource profile.
    //
    // The diagnostic counts read a mix of pre-QC (the raw beagle output) and post-QC (qc.modi)
    // files, matching the original. Reproduced as stdout only, captured in this task's own log.
    //
    // The per-chr QC'd VCF (qc.modi.vcf.gz) is published as well as converted to PLINK: it's what
    // the original 2_Axiom_LowPass/axiom_vs_lowpass.sh's bcftools isec call actually needs — 03_qc's
    // Axiom-vs-LowPass concordance study consumes it directly, VCF format, not the PLINK conversion.
    tag   "${chr}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e'

    input:
    tuple val(chr), path(beagle_vcf), path(beagle_csi)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools filter -e 'INFO/DR2<0.8' ${beagle_vcf} -Oz -o DA_AFFYimp_${chr}.qc.vcf.gz
    bcftools index DA_AFFYimp_${chr}.qc.vcf.gz

    bcftools annotate --set-id '%CHROM:%POS' DA_AFFYimp_${chr}.qc.vcf.gz -Oz -o DA_AFFYimp_${chr}.qc.modi.vcf.gz
    bcftools index DA_AFFYimp_${chr}.qc.modi.vcf.gz

    echo 'bcftools query -l DA_AFFYimp_${chr}.qc.modi.vcf.gz | wc -l'
    bcftools query -l DA_AFFYimp_${chr}.qc.modi.vcf.gz | wc -l

    echo 'bcftools view -H ${beagle_vcf} | wc -l'
    bcftools view -H ${beagle_vcf} | wc -l

    echo 'bcftools view -H DA_AFFYimp_${chr}.qc.modi.vcf.gz | wc -l'
    bcftools view -H DA_AFFYimp_${chr}.qc.modi.vcf.gz | wc -l

    plink --dog --vcf DA_AFFYimp_${chr}.qc.modi.vcf.gz --keep-allele-order --double-id --out DA_AFFYimp_${chr} --make-bed
    """

    stub:
    """
    touch DA_AFFYimp_${chr}.bed
    touch DA_AFFYimp_${chr}.bim
    touch DA_AFFYimp_${chr}.fam
    touch DA_AFFYimp_${chr}.qc.modi.vcf.gz
    touch DA_AFFYimp_${chr}.qc.modi.vcf.gz.csi
    """

    output:
    tuple val(chr), path("DA_AFFYimp_${chr}.bed"), path("DA_AFFYimp_${chr}.bim"), path("DA_AFFYimp_${chr}.fam"), emit: plink
    tuple val(chr), path("DA_AFFYimp_${chr}.qc.modi.vcf.gz"), path("DA_AFFYimp_${chr}.qc.modi.vcf.gz.csi"), emit: vcf
}


process MERGE_CHR_PLINK_AXIOM {
    // Cohort-level: merge all per-chr PLINK files into the final imputed Axiom dataset
    // Source: 02_data_processing/2_Axiom_imputation/03_axiomimp_to_plink.sh
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/axiom_liftover_impute:27db68ae33d39b0e'

    input:
    path(plink_files, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Build merge list (chr2-38; chr1 is the --bfile base)
    > merge_list.txt
    for i in \$(seq 2 38); do
        echo "DA_AFFYimp_chr\${i}" >> merge_list.txt
    done

    plink --allow-no-sex \\
        --dog \\
        --geno 0.05 \\
        --keep-allele-order \\
        --bfile DA_AFFYimp_chr1 \\
        --merge-list merge_list.txt \\
        --make-bed \\
        --out DA_AFFYimp_ALLCHR
    """

    stub:
    """
    touch DA_AFFYimp_ALLCHR.bed
    touch DA_AFFYimp_ALLCHR.bim
    touch DA_AFFYimp_ALLCHR.fam
    touch DA_AFFYimp_ALLCHR.log
    """

    output:
    path("DA_AFFYimp_ALLCHR.bed"), emit: bed
    path("DA_AFFYimp_ALLCHR.bim"), emit: bim
    path("DA_AFFYimp_ALLCHR.fam"), emit: fam
    path("DA_AFFYimp_ALLCHR.log"), emit: log
    tuple val("${task.process}"), val('plink'), eval("plink --version | sed 's/PLINK v//;s/ .*//'"), emit: versions, topic: versions
}
