// ============================================================================
// process_definitions.nf — 03. QC processes
//
// Generic VCF-comparison primitives (FILTER_PASS_VCF ... GENOTYPE_CONCORDANCE) are shared by
// all three concordance studies via the CONCORDANCE_ANALYSIS subworkflow in 03_qc.nf.
// DOWNSAMPLE_BAM is specific to the downsampling study.
//
// The *_DOWN processes below are forks of the equivalent 01_mapping processes with one change:
// a `meta` value (the downsample fraction) is threaded through the input/output tuples
// alongside chr/sample_id, since the unmodified 01_mapping processes only echo back chr or
// sample_id. Forking here rather than editing 01_mapping/process_definitions.nf directly keeps
// that already-reviewed pipeline code untouched, at the cost of script/stub bodies duplicated
// from there that need to stay in sync if the source processes change.
// ============================================================================

process FILTER_PASS_VCF {
    // Restrict a jointly-called VCF to PASS-filter sites
    // Source: 03_qc/1_High_LowPass/high_vs_lowpass.sh
    //
    // Resolved issue: the original ran `bcftools filter -i 'PASS'` immediately followed
    // by `bcftools view -f 'PASS'` into the same output filename. `-i/--include` takes a
    // boolean expression (documented and real-world usage always writes `FILTER="PASS"`,
    // never a bare literal), so `-i 'PASS'` is very likely invalid/non-functional syntax;
    // `-f/--apply-filters PASS` is the correct, standard way to keep only PASS sites and
    // is almost certainly what actually produced the archived output. Implemented here
    // using only the working `-f` form.
    tag   "${label}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(label), path(vcf), path(vcf_idx)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view -f 'PASS' ${vcf} -Oz -o ${label}.pass.vcf.gz
    bcftools index ${label}.pass.vcf.gz
    """

    stub:
    """
    touch ${label}.pass.vcf.gz
    touch ${label}.pass.vcf.gz.csi
    """

    output:
    tuple val(label), path("${label}.pass.vcf.gz"), path("${label}.pass.vcf.gz.csi"),                    emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process EXTRACT_SAMPLES_VCF {
    // Restrict a VCF to a fixed list of samples
    // Source: 03_qc/1_High_LowPass/high_vs_lowpass.sh; 03_qc/3_Downsampling/3_genotype_concordance.sh
    tag   "${label}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(label), path(vcf), path(vcf_idx), path(samples_list)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view -S ${samples_list} ${vcf} -Oz -o ${label}.extracted.vcf.gz
    bcftools index ${label}.extracted.vcf.gz
    """

    stub:
    """
    touch ${label}.extracted.vcf.gz
    touch ${label}.extracted.vcf.gz.csi
    """

    output:
    tuple val(label), path("${label}.extracted.vcf.gz"), path("${label}.extracted.vcf.gz.csi"),          emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process MAF_FILTER_VCF {
    // Exclude sites with MAF < 0.01
    // Source: 03_qc/1_High_LowPass/high_vs_lowpass.sh; 03_qc/2_Axiom_LowPass/axiom_vs_lowpass.sh
    tag   "${label}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(label), path(vcf), path(vcf_idx)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view -e 'MAF<0.01' ${vcf} -Oz -o ${label}.maf0.01.vcf.gz
    bcftools index ${label}.maf0.01.vcf.gz
    """

    stub:
    """
    touch ${label}.maf0.01.vcf.gz
    touch ${label}.maf0.01.vcf.gz.csi
    """

    output:
    tuple val(label), path("${label}.maf0.01.vcf.gz"), path("${label}.maf0.01.vcf.gz.csi"),              emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process CONCAT_CHR_VCFS {
    // Concatenate per-chromosome VCFs into one all-chromosome VCF
    // Source: 03_qc/2_Axiom_LowPass/axiom_vs_lowpass.sh
    // Chromosome order is derived with `sort -V` (natural/version sort: chr1, chr2, ..., chr10, ...)
    // rather than a fixed file list, so this works for any per-chr VCF set with parseable names.
    tag   "${label}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(label), path(chr_vcfs, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ls *.vcf.gz | sort -V > filelist.txt
    bcftools concat --file-list filelist.txt -Oz -o ${label}.ALLCHR.vcf.gz
    bcftools index ${label}.ALLCHR.vcf.gz
    """

    stub:
    """
    touch ${label}.ALLCHR.vcf.gz
    touch ${label}.ALLCHR.vcf.gz.csi
    """

    output:
    tuple val(label), path("${label}.ALLCHR.vcf.gz"), path("${label}.ALLCHR.vcf.gz.csi"),                emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process ISEC_COMMON_SITES {
    // Find sites present in both VCF A and VCF B
    // Source: 03_qc/1_High_LowPass/high_vs_lowpass.sh; 03_qc/2_Axiom_LowPass/axiom_vs_lowpass.sh
    tag   "${label}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(label), path(vcf_a), path(vcf_a_idx), path(vcf_b), path(vcf_b_idx)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools isec -p isec_out -n=2 -w1 ${vcf_a} ${vcf_b}
    cp isec_out/sites.txt ${label}.sites.txt
    """

    stub:
    """
    touch ${label}.sites.txt
    """

    output:
    tuple val(label), path("${label}.sites.txt"),                                                        emit: sites
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process RESTRICT_TO_SITES_VCF {
    // Restrict a VCF to a fixed site list and sort
    // Source: 03_qc/1_High_LowPass/high_vs_lowpass.sh; 03_qc/2_Axiom_LowPass/axiom_vs_lowpass.sh
    tag   "${label}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(label), path(vcf), path(vcf_idx), path(sites)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view -T ${sites} ${vcf} | bcftools sort -Oz -o ${label}.common.vcf.gz
    bcftools index ${label}.common.vcf.gz
    """

    stub:
    """
    touch ${label}.common.vcf.gz
    touch ${label}.common.vcf.gz.csi
    """

    output:
    tuple val(label), path("${label}.common.vcf.gz"), path("${label}.common.vcf.gz.csi"),                emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process UPDATE_VCF_DICT {
    // Reindex and rewrite a VCF's sequence dictionary to match the reference
    // Source: 03_qc/1_High_LowPass/high_vs_lowpass.sh; 03_qc/2_Axiom_LowPass/axiom_vs_lowpass.sh
    // Needed because the two VCFs entering GENOTYPE_CONCORDANCE were called with
    // differing contig sets in the original pipelines.
    tag   "${label}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(label), path(vcf), path(vcf_idx)
    path(ref_dict)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    gatk IndexFeatureFile -I ${vcf}
    gatk UpdateVCFSequenceDictionary \\
        -V ${vcf} \\
        --replace true \\
        --source-dictionary ${ref_dict} \\
        --output newcontig_${label}.vcf.gz
    bcftools index newcontig_${label}.vcf.gz
    """

    stub:
    """
    touch newcontig_${label}.vcf.gz
    touch newcontig_${label}.vcf.gz.csi
    """

    output:
    tuple val(label), path("newcontig_${label}.vcf.gz"), path("newcontig_${label}.vcf.gz.csi"),          emit: vcf
    tuple val("${task.process}"), val('gatk'),     eval("gatk --version | sed '/GATK/!d;s/.* v//'"),         emit: versions_gatk,     topic: versions
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process GENOTYPE_CONCORDANCE {
    // Compare one call/truth sample pair
    // Source: 03_qc/1_High_LowPass/high_vs_lowpass.sh; 03_qc/2_Axiom_LowPass/axiom_vs_lowpass.sh;
    //         03_qc/3_Downsampling/3_genotype_concordance.sh
    tag   "${pair_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89'

    input:
    tuple val(pair_id), val(call_sample), val(truth_sample),
          path(call_vcf), path(call_vcf_idx), path(truth_vcf), path(truth_vcf_idx)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    picard GenotypeConcordance \\
        -CALL_VCF ${call_vcf} \\
        -CALL_SAMPLE ${call_sample} \\
        -O ${pair_id}_concordance.vcf \\
        -TRUTH_VCF ${truth_vcf} \\
        -TRUTH_SAMPLE ${truth_sample}
    """

    stub:
    """
    touch ${pair_id}_concordance.vcf.genotype_concordance_summary_metrics
    touch ${pair_id}_concordance.vcf.genotype_concordance_contingency_metrics
    touch ${pair_id}_concordance.vcf.genotype_concordance_detail_metrics
    """

    output:
    tuple val(pair_id), path("${pair_id}_concordance.vcf.genotype_concordance_summary_metrics"),         emit: summary_metrics
    tuple val(pair_id), path("${pair_id}_concordance.vcf.genotype_concordance_contingency_metrics"),     emit: contingency_metrics
    tuple val(pair_id), path("${pair_id}_concordance.vcf.genotype_concordance_detail_metrics"),          emit: detail_metrics
    tuple val("${task.process}"), val('picard'), eval("picard MarkDuplicates --version 2>&1 | sed '/Version/!d; s/.*://'"), emit: versions_picard, topic: versions
}


process DOWNSAMPLE_BAM {
    // Downsample one HighPass BAM to a target fraction
    // Source: 03_qc/3_Downsampling/1_create_downsampled_bams.sh
    //
    // The downsampling probability P is computed here from a per-sample constant
    // (p_at_1x, taken verbatim from the original script's hardcoded per-sample table —
    // the depth measurement that originally produced it is not itself in the repo) rather
    // than hardcoding one P value per (sample, fraction) pair: P(fraction) = fraction * p_at_1x.
    // Verified against the original table, e.g. SRR13340532: p_at_1x=0.02443, and the
    // original's hardcoded 0.9x value 0.02198 = 0.02443 * 0.9 (rounded to 5 d.p.).
    tag   "${sample_id}:${fraction}"
    label 'process_low'
    container 'community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89'

    input:
    tuple val(sample_id), val(fraction), path(bam), path(bai), val(p_at_1x)

    when:
    task.ext.when == null || task.ext.when

    script:
    def target = fraction.replace('x', '') as Double
    def p      = target * (p_at_1x as Double)
    """
    picard DownsampleSam \\
        I=${bam} \\
        O=${sample_id}_${fraction}_ds.sorted.bam \\
        P=${p}
    samtools index ${sample_id}_${fraction}_ds.sorted.bam
    """

    stub:
    """
    touch ${sample_id}_${fraction}_ds.sorted.bam
    touch ${sample_id}_${fraction}_ds.sorted.bam.bai
    """

    output:
    tuple val(sample_id), val(fraction),
          path("${sample_id}_${fraction}_ds.sorted.bam"),
          path("${sample_id}_${fraction}_ds.sorted.bam.bai"),                                            emit: bam
    tuple val("${task.process}"), val('picard'),   eval("picard MarkDuplicates --version 2>&1 | sed '/Version/!d; s/.*://'"), emit: versions_picard,   topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/.* //'"),                          emit: versions_samtools, topic: versions
}


process MERGE_BQSR_GVCF_DOWN {
    // Fork of 01_mapping's DOG10K_MERGE_MARKDUPS_BQSR_GVCF — adds a `meta` (downsample
    // fraction) passthrough. See file header for why this is a fork, not a shared include.
    // Per-sample: merge all runs -> MarkDups -> BQSR -> HaplotypeCaller -> gVCF
    tag   "${meta}:${sample_id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89'

    input:
    tuple val(meta), val(sample_id), path(bams, arity: '1..*'), path(bais, arity: '1..*')
    tuple path(assembly_ref), path(assembly_ref_idx, arity: '1..*')
    tuple path(known_variants), path(known_variants_tbi, arity: '1')
    path(chunks_dir)

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${sample_id}"
    """
    ## merge all runs for this sample
    ## sorted BAMs are staged from DOG10K_MAPPING and named \${sample_id}_\${run_id}.sorted.bam
    ls ${prefix}_*.sorted.bam > ${prefix}.bamlist
    if [ \$(cat ${prefix}.bamlist | wc -l) -gt 1 ]; then
        samtools merge -@ ${task.cpus} -b ${prefix}.bamlist ${prefix}.sorted.merged.bam
    else
        mv \$(cat ${prefix}.bamlist) ${prefix}.sorted.merged.bam
    fi
    samtools index ${prefix}.sorted.merged.bam

    ## extract the unmapped reads
    samtools view -f 4 ${prefix}.sorted.merged.bam \\
        | perl -lane '\$line=\$_;if(\$F[0]=~/@/){print}else{print if \$F[2] eq "*"}' \\
        | samtools view -b > ${prefix}.unmapped.bam

    ## split the bam files, and generate the BQSR table for each chunk
    chunk_process_s1() {
        local chunk_bed=\$1
        local chunk_i=\$(basename \$chunk_bed ".bed")

        ## split the bam
        samtools view -hb -L \$chunk_bed ${prefix}.sorted.merged.bam \\
            > ${prefix}.sorted.merged.\${chunk_i}.bam

        ## mark duplicates
        picard MarkDuplicates \\
            METRICS_FILE=${prefix}.sorted.merged.\${chunk_i}.matrix \\
            INPUT=${prefix}.sorted.merged.\${chunk_i}.bam \\
            OUTPUT=${prefix}.sorted.merged.MarkDups.\${chunk_i}.bam

        samtools index ${prefix}.sorted.merged.MarkDups.\${chunk_i}.bam

        ## BQSR
        # Generate the first pass BQSR table file
        gatk --java-options "-Xmx4G" BaseRecalibrator \\
            -R ${assembly_ref} \\
            -I ${prefix}.sorted.merged.MarkDups.\${chunk_i}.bam \\
            -L \$chunk_bed \\
            --known-sites ${known_variants} \\
            -O ${prefix}.sorted.merged.MarkDups.\${chunk_i}.table
    }

    ## apply the BQSR for each chunk
    chunk_process_s2() {
        local chunk_bed=\$1
        local chunk_i=\$(basename \$chunk_bed ".bed")

        # Apply BQSR
        gatk --java-options "-Xmx4G" ApplyBQSR \\
            -R ${assembly_ref} \\
            -I ${prefix}.sorted.merged.MarkDups.\${chunk_i}.bam \\
            -L \$chunk_bed \\
            -bqsr ${prefix}.BQSR.reports.table \\
            --preserve-qscores-less-than 6 \\
            --static-quantized-quals 10 \\
            --static-quantized-quals 20 \\
            --static-quantized-quals 30 \\
            -O ${prefix}.sorted.merged.MarkDups.BQSR.\${chunk_i}.bam
    }

    chunk_process_s3() {
        local chunk_bed=\$1
        local chunk_i=\$(basename \$chunk_bed ".bed")

        ## GVCF
        gatk --java-options "-Xmx4G" HaplotypeCaller \\
            -R ${assembly_ref} \\
            -ERC GVCF \\
            -L \$chunk_bed \\
            -OVI \\
            -I ${prefix}.sorted.merged.MarkDups.BQSR.bam \\
            -O ${prefix}.sorted.merged.MarkDups.BQSR.\${chunk_i}.g.vcf.gz

        tabix ${prefix}.sorted.merged.MarkDups.BQSR.\${chunk_i}.g.vcf.gz
    }

    ## split the gvcf for each chromosome, in order to merge them using GatherVCFs
    chunk_process_s4() {
        local chunk_num=\$1
        local chunk_i="chunk_\${chunk_num}"

        tabix --list-chroms ${prefix}.sorted.merged.MarkDups.BQSR.\${chunk_i}.g.vcf.gz \\
            > \${chunk_i}.list
        zcat ${prefix}.sorted.merged.MarkDups.BQSR.\${chunk_i}.g.vcf.gz \\
            | grep "^#" > \${chunk_i}.header

        while IFS= read -r chrom; do
            tabix ${prefix}.sorted.merged.MarkDups.BQSR.\${chunk_i}.g.vcf.gz \$chrom \\
                | cat \${chunk_i}.header - \\
                | bgzip > ${prefix}.sorted.merged.MarkDups.BQSR.sc.\${chrom}.g.vcf.gz
            tabix -p vcf ${prefix}.sorted.merged.MarkDups.BQSR.sc.\${chrom}.g.vcf.gz
        done < \${chunk_i}.list
    }

    ## split the genome into 10 chunks, and process s1 for each chunk
    for chunk_bed in ${chunks_dir}/*.bed; do
        chunk_process_s1 \$chunk_bed &
    done
    wait

    ## combine the BQSR reports from each chunk
    ls ${prefix}.sorted.merged.MarkDups.*.table > ${prefix}.BQSR.reports.list
    #
    gatk --java-options "-Xmx6G" GatherBQSRReports \\
        -I ${prefix}.BQSR.reports.list \\
        -O ${prefix}.BQSR.reports.table
    xargs -a ${prefix}.BQSR.reports.list rm -f
    #
    rm -f ${prefix}.sorted.merged.bam
    rm -f ${prefix}.sorted.merged.bam.bai
    rm -f ${prefix}.sorted.merged.chunk_*.matrix
    rm -f ${prefix}.sorted.merged.chunk_*.bam
    rm -f ${prefix}.sorted.merged.chunk_*.bam.bai
    #
    ## apply BQSR with combined table, and continue with s2 process
    for chunk_bed in ${chunks_dir}/*.bed; do
        chunk_process_s2 \$chunk_bed &
    done
    wait

    ## merge the bam chunks
    ls ${prefix}.sorted.merged.MarkDups.BQSR.chunk_*.bam > chunks.bamlist
    echo "${prefix}.unmapped.bam" >> chunks.bamlist
    samtools merge -@ ${task.cpus} -fb chunks.bamlist \\
        ${prefix}.sorted.merged.MarkDups.BQSR.bam
    samtools index ${prefix}.sorted.merged.MarkDups.BQSR.bam

    rm -f chunks.bamlist
    rm -f ${prefix}.sorted.merged.MarkDups.chunk_*.bam
    rm -f ${prefix}.sorted.merged.MarkDups.chunk_*.bam.bai
    rm -f ${prefix}.sorted.merged.MarkDups.BQSR.chunk_*.bam
    rm -f ${prefix}.sorted.merged.MarkDups.BQSR.chunk_*.bai
    rm -f ${prefix}.unmapped.bam

    ## generate the g.vcf file for each chunk
    for chunk_bed in ${chunks_dir}/*.bed; do
        chunk_process_s3 \$chunk_bed &
    done
    wait

    ## no need to sort the contigs in chunk_11, just split the first 10 chunks
    for chunk_num in {1..10}; do
        chunk_process_s4 \$chunk_num &
    done
    wait

    ## merge the g.vcf files
    gatk --java-options "-Xmx6G" GatherVcfs \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr1.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr2.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr3.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr4.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr5.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr6.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr7.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr8.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr9.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr10.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr11.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr12.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr13.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr14.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr15.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr16.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr17.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr18.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr19.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr20.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr21.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr22.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr23.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr24.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr25.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr26.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr27.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr28.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr29.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr30.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr31.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr32.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr33.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr34.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr35.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr36.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr37.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr38.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chrX.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chrY_NC_051844.1.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chrM.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chrY_unplaced_NW_024010443.1.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.sc.chrY_unplaced_NW_024010444.1.g.vcf.gz \\
        -I ${prefix}.sorted.merged.MarkDups.BQSR.chunk_11.g.vcf.gz \\
        -O ${prefix}.sorted.merged.MarkDups.BQSR.g.vcf.gz

    tabix -p vcf ${prefix}.sorted.merged.MarkDups.BQSR.g.vcf.gz

    ## clean the files
    rm -f chunk_*.list chunk_*.header
    rm -f ${prefix}.BQSR.reports.list
    rm -f ${prefix}.sorted.merged.MarkDups.BQSR.chunk_*.g.vcf.gz
    rm -f ${prefix}.sorted.merged.MarkDups.BQSR.chunk_*.g.vcf.gz.tbi
    rm -f ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr*.g.vcf.gz
    rm -f ${prefix}.sorted.merged.MarkDups.BQSR.sc.chr*.g.vcf.gz.tbi
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_id}"
    """
    touch ${prefix}.sorted.merged.MarkDups.BQSR.bam
    touch ${prefix}.sorted.merged.MarkDups.BQSR.bam.bai
    touch ${prefix}.sorted.merged.MarkDups.BQSR.g.vcf.gz
    touch ${prefix}.sorted.merged.MarkDups.BQSR.g.vcf.gz.tbi
    """

    output:
    tuple val(meta), val(sample_id),
          path("*.sorted.merged.MarkDups.BQSR.bam"),
          path("*.sorted.merged.MarkDups.BQSR.bam.bai"),                                                                      emit: bam
    tuple val(meta), val(sample_id),
          path("*.sorted.merged.MarkDups.BQSR.g.vcf.gz"),
          path("*.sorted.merged.MarkDups.BQSR.g.vcf.gz.tbi"),                                                                 emit: gvcf
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/.* //'"),                           emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('picard'),   eval("picard MarkDuplicates --version 2>&1 | sed '/Version/!d; s/.*://'"), emit: versions_picard,   topic: versions
    tuple val("${task.process}"), val('gatk'),     eval("gatk --version | sed '/GATK/!d;s/.* v//'"),                          emit: versions_gatk,     topic: versions
    tuple val("${task.process}"), val('tabix'),    eval("tabix --version | sed '/tabix/!d; s/.* //'"),                        emit: versions_tabix,    topic: versions
}


process GLIMPSE_LIKELIHOODS_DOWN {
    // Fork of 01_mapping's GLIMPSE_LIKELIHOODS — adds a `meta` (downsample fraction) passthrough.
    // Per-chr: compute genotype likelihoods from a bundle of remapped downsampled BAMs (bcftools only)
    tag   "${meta}:${chr}"
    label 'process_high'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(meta), val(chr), path(sites_vcf), path(sites_vcf_csi), path(sites_tsv), path(sites_tsv_tbi),
          path(bams), path(bais)
    tuple path(assembly_ref), path(assembly_ref_idx, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Create BAM list from staged files
    ls *.sorted.merged.MarkDups.BQSR.bam > bam_list.txt

    bcftools mpileup -f ${assembly_ref} -I -E -a 'FORMAT/DP' \\
        -T ${sites_vcf} \\
        -r ${chr} \\
        --bam-list bam_list.txt \\
        -Ou \\
        | bcftools call -Aim -C alleles \\
        -T ${sites_tsv} \\
        -Oz -o ${chr}.GL.DownSampled.${meta}.vcf.gz

    bcftools index -f ${chr}.GL.DownSampled.${meta}.vcf.gz
    """

    stub:
    """
    touch ${chr}.GL.DownSampled.${meta}.vcf.gz
    touch ${chr}.GL.DownSampled.${meta}.vcf.gz.csi
    """

    output:
    tuple val(meta), val(chr),
          path("${chr}.GL.DownSampled.${meta}.vcf.gz"),
          path("${chr}.GL.DownSampled.${meta}.vcf.gz.csi"),                                                  emit: gl_vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process GLIMPSE_IMPUTE_CHUNKS_DOWN {
    // Fork of 01_mapping's GLIMPSE_IMPUTE_CHUNKS — adds a `meta` (downsample fraction) passthrough.
    // Per-chunk: phase and impute one genomic chunk (GLIMPSE only)
    tag   "${meta}:${chr}:${chunk_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/glimpse-bio:1.1.1--7517c3cd783f15ab'

    input:
    tuple val(meta), val(chr), val(chunk_id), val(input_region), val(output_region),
          path(gl_vcf), path(gl_vcf_csi)
    tuple path(ref_panel), path(ref_panel_csi, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf -v PADDED_ID "%02d" ${chunk_id}

    GLIMPSE_phase \\
        --input ${gl_vcf} \\
        --reference ${ref_panel} \\
        --input-region ${input_region} \\
        --output-region ${output_region} \\
        --output Down_${meta}.${chr}.imputed.\${PADDED_ID}.bcf
    """

    stub:
    """
    printf -v PADDED_ID "%02d" ${chunk_id}
    touch Down_${meta}.${chr}.imputed.\${PADDED_ID}.bcf
    """

    output:
    tuple val(meta), val(chr), val(chunk_id),
          path("Down_${meta}.${chr}.imputed.*.bcf", arity: '1'),                                                emit: imputed_bcf
    tuple val("${task.process}"), val('GLIMPSE'),  eval("GLIMPSE_phase --help | sed '/Version/!d; s/.* : //'"), emit: versions_glimpse,  topic: versions
}


process INDEX_BCFS_DOWN {
    // Fork of 01_mapping's INDEX_BCFS — adds a `meta` (downsample fraction) passthrough.
    // Per-chunk: index imputed BCF files with bcftools (bcftools only)
    tag   "${meta}:${chr}:${chunk_id}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(meta), val(chr), val(chunk_id), path(bcf, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools index -f ${bcf}
    """

    stub:
    """
    touch ${bcf}.csi
    """

    output:
    tuple val(meta), val(chr), val(chunk_id),
          path(bcf, arity: '1'),
          path("${bcf}.csi", arity: '1'),                                                                    emit: imputed
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process GLIMPSE_LIGATE_DOWN {
    // Fork of 01_mapping's GLIMPSE_LIGATE — adds a `meta` (downsample fraction) passthrough.
    // Per-chr: ligate imputed chunks -> one merged BCF per chromosome (GLIMPSE only)
    tag   "${meta}:${chr}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/glimpse-bio:1.1.1--7517c3cd783f15ab'

    input:
    tuple val(meta), val(chr), path(chunk_bcfs, arity: '1..*'), path(chunk_csis, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Sort chunks by filename for correct ligation order.
    ls Down_${meta}.${chr}.imputed.*.bcf | sort > list.${chr}.txt

    GLIMPSE_ligate \\
        --input list.${chr}.txt \\
        --output Down_${meta}.${chr}.merged.bcf
    """

    stub:
    """
    touch Down_${meta}.${chr}.merged.bcf
    """

    output:
    tuple val(meta), val(chr),
          path("Down_${meta}.${chr}.merged.bcf"),                                                                emit: merged_bcf
    tuple val("${task.process}"), val('GLIMPSE'), eval("GLIMPSE_ligate --help | sed '/Version/!d; s/.* : //'"), emit: versions_glimpse,  topic: versions
}


process BCF_VCF_DOWN {
    // Fork of 01_mapping's BCF_VCF — adds a `meta` (downsample fraction) passthrough.
    // Per-chr: convert merged BCF -> bgzipped VCF + index (bcftools only)
    tag   "${meta}:${chr}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(meta), val(chr), path(bcfs, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view ${bcfs} -Oz -o Down_${meta}.${chr}.merged.vcf.gz
    bcftools index -f Down_${meta}.${chr}.merged.vcf.gz
    """

    stub:
    """
    touch Down_${meta}.${chr}.merged.vcf.gz
    touch Down_${meta}.${chr}.merged.vcf.gz.csi
    """

    output:
    tuple val(meta), val(chr),
          path("Down_${meta}.${chr}.merged.vcf.gz"),
          path("Down_${meta}.${chr}.merged.vcf.gz.csi"),                                                     emit: merged
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process GLIMPSE_QC_FILTER_DOWN {
    // Fork of 01_mapping's GLIMPSE_QC_FILTER — adds a `meta` (downsample fraction) passthrough.
    // Per-chr: filter on imputation INFO score and annotate variant IDs (bcftools only)
    tag   "${meta}:${chr}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(meta), val(chr), path(merged_vcf), path(merged_vcf_idx)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools filter -e 'INFO/INFO<0.8' ${merged_vcf} -Oz -o Down_${meta}.${chr}.merged.qc.vcf.gz
    bcftools annotate --set-id '%CHROM:%POS' Down_${meta}.${chr}.merged.qc.vcf.gz -Oz -o Down_${meta}.${chr}.merged.qc_modi.vcf.gz
    bcftools index -f Down_${meta}.${chr}.merged.qc_modi.vcf.gz
    """

    stub:
    """
    touch Down_${meta}.${chr}.merged.qc_modi.vcf.gz
    touch Down_${meta}.${chr}.merged.qc_modi.vcf.gz.csi
    """

    output:
    tuple val(meta), val(chr),
          path("Down_${meta}.${chr}.merged.qc_modi.vcf.gz"),
          path("Down_${meta}.${chr}.merged.qc_modi.vcf.gz.csi"),                                             emit: qc_vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}
