// ============================================================================
// process_definitions.nf — 01_mapping / 2_cohort_genotyping
//
// Cohort-wide processes: joint SNP calling (HighPass) and GLIMPSE imputation
// (LowPass). This pipeline is invoked exactly ONCE over the complete sample set,
// after every 1_fetch_map batch (both tracks) has completed — see 2_cohort_
// genotyping.nf and archive/conversion_notes/01_mapping.md. Running it per-batch
// instead would silently change the joint-calling / imputation output, not just
// be a performance difference.
//
// BWA_MEM2_INDEX is duplicated from 1_fetch_map/process_definitions.nf (same
// process, same container) so this pipeline can build its own ch_assembly_ref
// value channel exactly as the original single pipeline did — see conversion
// notes for the redundant-reindexing tradeoff this implies.
//
// All other process script bodies below are unchanged from the original
// code/01_mapping/process_definitions.nf — only which pipeline file they live
// in changed.
// ============================================================================

process BWA_MEM2_INDEX {
    tag   "reference"
    label 'process_high'
    container 'community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89'

    input:
    path(fasta, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bwa-mem2 index ${fasta}
    """

    stub:
    """
    touch ${fasta}.0123
    touch ${fasta}.amb
    touch ${fasta}.ann
    touch ${fasta}.bwt.2bit.64
    touch ${fasta}.pac
    """

    output:
    tuple path(fasta, arity: '1'),
          path("${fasta}.*", arity: '1..*'),                                                              emit: index
    tuple val("${task.process}"), val('bwa-mem2'), eval("bwa-mem2 version"),                              emit: versions_bwamem2, topic: versions
}




process EXTRACT_DEPTH {
    // Cohort-level: aggregate per-sample depth files → combined TSV
    // Source: 01_mapping/1_High-LowPass_mapping/extract_depth_autosomes_and_x.sh
    tag   "cohort"
    label 'process_single'
    container 'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156'

    input:
    path(depth_files, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Output file that will contain the combined results
    output_file="combined_knownsites.depth.txt"

    # Clear the output file if it already exists
    > "\$output_file"

    # Loop over all matching files, regardless of whether the sample ID starts with
    # letters, numbers, or a mix of both
    for file in *.sorted.merged.MarkDups.BQSR.bam.knownsites.depth.txt; do
      # Skip the loop if no files match the pattern
      [ -e "\$file" ] || continue

      # Remove the fixed suffix from the filename to get the sample ID
      sample_id=\$(basename "\$file" | sed 's/\\.sorted\\.merged\\.MarkDups\\.BQSR\\.bam\\.knownsites\\.depth\\.txt\$//')

      # Extract the 3rd column from rows 2 and 3, then join them with tabs
      values=\$(awk 'NR==2 || NR==3 {print \$3}' "\$file" | paste -sd '\\t' -)

      # Write sample ID and values to the output file
      printf '%s\\t%s\\n' "\$sample_id" "\$values" >> "\$output_file"
    done
    """

    stub:
    """
    touch combined_knownsites.depth.txt
    """

    output:
    path("combined_knownsites.depth.txt"), emit: depth_stats
}


process GENOTYPE_CALLING {
    // Cohort-level: joint calling + SNP hard filtering across all HighPass samples
    // Source: 01_mapping/2_HighPass_SNPfiltering/genomicsDB_joincall_SNPIDandHardFilter.sh
    //
    // [TODO] VariantRecalibrator section is commented out in the original and here:
    //   (1) Original references raw_output.vcf.gz but MergeVcfs outputs raw_output2.vcf.gz
    //   (2) ApplyVQSR was already commented out in the original
    //   (3) Recalibration resource VCF not wired as an input
    //   The pipeline proceeds directly to hard filtering, as in the original.
    //
    // [TODO] Fixed: original `bcftools annotate` had `Oz` (positional) not `-Oz` (flag)
    tag   "cohort"
    label 'process_high'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    path(gvcfs, arity: '1..*')
    path(gvcf_tbis, arity: '1..*')
    tuple path(assembly_ref), path(assembly_ref_idx, arity: '1..*')
    path(chunks_dir)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Create cohort.sample_map from staged gVCFs
    ## uses cohort.sample_map — put this in the run folder
    > cohort.sample_map
    for gvcf in *.sorted.merged.MarkDups.BQSR.g.vcf.gz; do
        id=\$(basename "\$gvcf" .sorted.merged.MarkDups.BQSR.g.vcf.gz)
        printf '%s\\t%s\\n' "\$id" "\$gvcf" >> cohort.sample_map
    done

    chunk_process_s1() {
        local i=\$1
        local chunk_i=\$(basename \$i ".bed")

        gatk --java-options "-Xmx4G" GenomicsDBImport \\
            --genomicsdb-workspace-path "my_database.\${chunk_i}" \\
            --sample-name-map cohort.sample_map \\
            --batch-size 50 \\
            --genomicsdb-shared-posixfs-optimizations true \\
            -L \$i

        gatk --java-options "-Xmx4G" GenotypeGVCFs \\
            -R ${assembly_ref} \\
            -V "gendb://my_database.\${chunk_i}" \\
            -L \$i \\
            -O "\${chunk_i}.raw_output.vcf.gz"
    }

    chunk_process_s2() {
        local chunk_i="chunk_\$1"
        tabix --list-chroms "\${chunk_i}.raw_output.vcf.gz" > "\${chunk_i}.list"
        zcat "\${chunk_i}.raw_output.vcf.gz" | grep "^#" > "\${chunk_i}.header"
        while IFS= read -r line; do
            tabix "\${chunk_i}.raw_output.vcf.gz" \$line \\
                | cat "\${chunk_i}.header" - \\
                | bgzip > "\${line}.raw_output.vcf.gz"
            tabix -p vcf "\${line}.raw_output.vcf.gz"
        done < "\${chunk_i}.list"
    }

    for j in ${chunks_dir}/*.bed; do
        (
        chunk_process_s1 \$j
        )&
    done
    wait

    for j in {1..10}; do
        (
        chunk_process_s2 \$j
        )&
    done
    wait

    gatk --java-options "-Xmx6G" MergeVcfs \\
        -I chr1.raw_output.vcf.gz \\
        -I chr2.raw_output.vcf.gz \\
        -I chr3.raw_output.vcf.gz \\
        -I chr4.raw_output.vcf.gz \\
        -I chr5.raw_output.vcf.gz \\
        -I chr6.raw_output.vcf.gz \\
        -I chr7.raw_output.vcf.gz \\
        -I chr8.raw_output.vcf.gz \\
        -I chr9.raw_output.vcf.gz \\
        -I chr10.raw_output.vcf.gz \\
        -I chr11.raw_output.vcf.gz \\
        -I chr12.raw_output.vcf.gz \\
        -I chr13.raw_output.vcf.gz \\
        -I chr14.raw_output.vcf.gz \\
        -I chr15.raw_output.vcf.gz \\
        -I chr16.raw_output.vcf.gz \\
        -I chr17.raw_output.vcf.gz \\
        -I chr18.raw_output.vcf.gz \\
        -I chr19.raw_output.vcf.gz \\
        -I chr20.raw_output.vcf.gz \\
        -I chr21.raw_output.vcf.gz \\
        -I chr22.raw_output.vcf.gz \\
        -I chr23.raw_output.vcf.gz \\
        -I chr24.raw_output.vcf.gz \\
        -I chr25.raw_output.vcf.gz \\
        -I chr26.raw_output.vcf.gz \\
        -I chr27.raw_output.vcf.gz \\
        -I chr28.raw_output.vcf.gz \\
        -I chr29.raw_output.vcf.gz \\
        -I chr30.raw_output.vcf.gz \\
        -I chr31.raw_output.vcf.gz \\
        -I chr32.raw_output.vcf.gz \\
        -I chr33.raw_output.vcf.gz \\
        -I chr34.raw_output.vcf.gz \\
        -I chr35.raw_output.vcf.gz \\
        -I chr36.raw_output.vcf.gz \\
        -I chr37.raw_output.vcf.gz \\
        -I chr38.raw_output.vcf.gz \\
        -I chrX.raw_output.vcf.gz \\
        -I chrY_NC_051844.1.raw_output.vcf.gz \\
        -I chrM.raw_output.vcf.gz \\
        -I chrY_unplaced_NW_024010443.1.raw_output.vcf.gz \\
        -I chrY_unplaced_NW_024010444.1.raw_output.vcf.gz \\
        -I chunk_11.raw_output.vcf.gz \\
        -O raw_output2.vcf.gz

    #gatk --java-options "-Xmx5g -Xms5g" VariantRecalibrator \
    #-R ${assembly_ref} \
    #-V raw_output2.vcf.gz \
    #--resource:array,known=false,training=true,truth=true,prior=12.0 [FILL ME IN: recalibration_resource.vcf.gz] \
    #-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
    #-mode SNP \
    #-O snp.output.recal \
    #--tranches-file snp.output.tranches \
    #--rscript-file snp.output.plots.R \
    #--trust-all-polymorphic true

    #gatk --java-options "-Xmx5g -Xms5g" ApplyVQSR \
    #-V raw_output2.vcf.gz \
    #--recal-file snp.output.recal \
    #--tranches-file snp.output.tranches \
    #--truth-sensitivity-filter-level 99.0 \
    #--create-output-variant-index true \
    #-mode SNP \
    #-O snp.recal99.vcf.gz

    ####SNPIDandFilterSteps:
    # STEP 1 ## Annotate SNPs
    ## [TODO] Fixed: original had positional `Oz` instead of flag `-Oz`
    bcftools annotate --set-id '%CHROM\\:%POS' raw_output2.vcf.gz -Oz -o raw_output.ann.id.vcf.gz
    tabix raw_output.ann.id.vcf.gz

    ## STEP 4 ## Hard filters
    # Run HF SNP on All - as per GATK
    gatk SelectVariants -V raw_output.ann.id.vcf.gz -select-type SNP -O SNP.raw_output.ann.id.vcf.gz

    gatk VariantFiltration -V SNP.raw_output.ann.id.vcf.gz \\
        -filter "QD < 2.0" --filter-name "QD2" \\
        -filter "QUAL < 30.0" --filter-name "QUAL30" \\
        -filter "SOR > 3.0" --filter-name "SOR3" \\
        -filter "FS > 60.0" --filter-name "FS60" \\
        -filter "MQ < 40.0" --filter-name "MQ40" \\
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \\
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \\
        -O SNP.HF.ann.id.vcf.gz
    tabix SNP.HF.ann.id.vcf.gz

    #zcat SNP.HF.ann.id.vcf.gz| grep -v "#" | less -S

    ## cleanup intermediates
    rm -f chunk_*.list chunk_*.header
    rm -f *.raw_output.vcf.gz *.raw_output.vcf.gz.tbi
    rm -f raw_output2.vcf.gz raw_output2.vcf.gz.tbi
    rm -f raw_output.ann.id.vcf.gz raw_output.ann.id.vcf.gz.tbi
    rm -f SNP.raw_output.ann.id.vcf.gz SNP.raw_output.ann.id.vcf.gz.tbi SNP.raw_output.ann.id.vcf.gz.idx
    rm -rf my_database.*
    """

    stub:
    """
    touch SNP.HF.ann.id.vcf.gz
    touch SNP.HF.ann.id.vcf.gz.tbi
    """

    output:
    path("SNP.HF.ann.id.vcf.gz"),                                                                            emit: vcf
    path("SNP.HF.ann.id.vcf.gz.tbi"),                                                                        emit: vcf_tbi
    tuple val("${task.process}"), val('gatk'),     eval("gatk --version | sed '/GATK/!d;s/.* v//'"),         emit: versions_gatk,     topic: versions
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
    tuple val("${task.process}"), val('tabix'),    eval("tabix --version | sed '/tabix/!d; s/.* //'"),       emit: versions_tabix,    topic: versions
}


process GLIMPSE_DEFINE_CHUNKS {
    // Per-chr, one-time: define imputation chunks from phased reference panel
    // Source: 01_mapping/3_LowPass_imputation/scripts/DoOnce_GlimpseChunks.sh
    tag   "${chr}"
    label 'process_single'
    container 'community.wave.seqera.io/library/glimpse-bio:1.1.1--7517c3cd783f15ab'

    input:
    val(chr)
    tuple path(ref_panel), path(ref_panel_csi, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    GLIMPSE_chunk \\
        --input ${ref_panel} \\
        --region ${chr} \\
        --output chunks_${chr}.txt
    """

    stub:
    """
    echo "0 ${chr} ${chr}:1-2000000 ${chr}:1-1000000" > chunks_${chr}.txt
    """

    output:
    tuple val(chr), path("chunks_${chr}.txt", arity: '1'),                                                     emit: chunks
    tuple val("${task.process}"), val('GLIMPSE'), eval("GLIMPSE_chunk --help | sed '/Version/!d; s/.* : //'"), emit: versions, topic: versions
}


process GLIMPSE_VARIABLE_SITES {
    // Per-chr: extract polymorphic biallelic SNPs from phased Dog10K panel (bcftools + tabix only)
    // Source: 01_mapping/3_LowPass_imputation/scripts/01_VariableSites.sh
    // Note: bcftools-only process — uses highpass_snpfiltering container
    tag   "${chr}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    val(chr)
    tuple path(ref_panel), path(ref_panel_csi, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # We need a VCF/BCF file containing only sites to tell BCFtools at which positions
    # to make a genotype call. Since BCFtools does not compute correctly genotype
    # likelihoods for indels, here we only focus on SNPs (however, GLIMPSE can impute
    # any type of variants as soon as it is bi-allelic and has GLs being properly defined).
    # To perform the extraction from one chromosome of a reference panel, run BCFtools as follows:

    # -G -m 2 -M 2 -v snps:
    #   -G  group for HWE filter
    #   -m  multiallelic caller for rare variant calling
    #   -M  keep masked
    #   -v snps  keep variant sites only
    bcftools view -G -m 2 -M 2 -v snps -r ${chr} ${ref_panel} \\
        -Oz -o ${chr}.Dog10K.Phased.snp.sites.vcf.gz
    bcftools index -f ${chr}.Dog10K.Phased.snp.sites.vcf.gz

    # Then, convert the output file into TSV format and index using tabix (requires htslib in PATH):
    bcftools query -f'%CHROM\\t%POS\\t%REF,%ALT\\n' ${chr}.Dog10K.Phased.snp.sites.vcf.gz \\
        | bgzip -c > ${chr}.Dog10K.Phased.snp.sites.tsv.gz
    tabix -s1 -b2 -e2 ${chr}.Dog10K.Phased.snp.sites.tsv.gz
    """

    stub:
    """
    touch ${chr}.Dog10K.Phased.snp.sites.vcf.gz
    touch ${chr}.Dog10K.Phased.snp.sites.vcf.gz.csi
    touch ${chr}.Dog10K.Phased.snp.sites.tsv.gz
    touch ${chr}.Dog10K.Phased.snp.sites.tsv.gz.tbi
    """

    output:
    tuple val(chr),
          path("${chr}.Dog10K.Phased.snp.sites.vcf.gz"),
          path("${chr}.Dog10K.Phased.snp.sites.vcf.gz.csi"),
          path("${chr}.Dog10K.Phased.snp.sites.tsv.gz"),
          path("${chr}.Dog10K.Phased.snp.sites.tsv.gz.tbi"),                                                 emit: sites
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
    tuple val("${task.process}"), val('tabix'),    eval("tabix --version | sed '/tabix/!d; s/.* //'"),       emit: versions_tabix,    topic: versions
}


process GLIMPSE_LIKELIHOODS {
    // Per-chr, cohort-level: compute genotype likelihoods from all LowPass BAMs (bcftools only)
    // Source: 01_mapping/3_LowPass_imputation/scripts/02_GenotypeLiklihood.sh
    // Note! This is the longest step — can take ~2 days for chr1 on ~1800 samples
    // Note: bcftools-only process — uses highpass_snpfiltering container
    tag   "${chr}"
    label 'process_high'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(chr), path(sites_vcf), path(sites_vcf_csi), path(sites_tsv), path(sites_tsv_tbi),
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
        -Oz -o ${chr}.GL.DA_LowPass.vcf.gz

    bcftools index -f ${chr}.GL.DA_LowPass.vcf.gz
    """

    stub:
    """
    touch ${chr}.GL.DA_LowPass.vcf.gz
    touch ${chr}.GL.DA_LowPass.vcf.gz.csi
    """

    output:
    tuple val(chr),
          path("${chr}.GL.DA_LowPass.vcf.gz"),
          path("${chr}.GL.DA_LowPass.vcf.gz.csi"),                                                           emit: gl_vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process GLIMPSE_IMPUTE_CHUNKS {
    // Per-chunk: phase and impute one genomic chunk (GLIMPSE only)
    // Source: 01_mapping/3_LowPass_imputation/scripts/03_ImputeAndPhaseChrsByChunks.sh
    //
    // The input for the script is one line with the chunk coordinates defined.
    // Chunk IDs are zero-padded to match the original printf -v ID "%02d" pattern.
    // [SPLIT] bcftools index moved to INDEX_BCFS process
    tag   "${chr}:${chunk_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/glimpse-bio:1.1.1--7517c3cd783f15ab'

    input:
    tuple val(chr), val(chunk_id), val(input_region), val(output_region),
          path(gl_vcf), path(gl_vcf_csi)
    tuple path(ref_panel), path(ref_panel_csi, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    printf -v PADDED_ID "%02d" ${chunk_id}

    #OBS!! sample naming matches DA.LowPass.vcf.gz convention from original
    GLIMPSE_phase \\
        --input ${gl_vcf} \\
        --reference ${ref_panel} \\
        --input-region ${input_region} \\
        --output-region ${output_region} \\
        --output DA.${chr}.imputed.\${PADDED_ID}.bcf
    """

    stub:
    """
    printf -v PADDED_ID "%02d" ${chunk_id}
    touch DA.${chr}.imputed.\${PADDED_ID}.bcf
    """

    output:
    tuple val(chr), val(chunk_id),
          path("DA.${chr}.imputed.*.bcf", arity: '1'),                                                          emit: imputed_bcf
    tuple val("${task.process}"), val('GLIMPSE'),  eval("GLIMPSE_phase --help | sed '/Version/!d; s/.* : //'"), emit: versions_glimpse,  topic: versions
}


process INDEX_BCFS {
    // Per-chunk: index imputed BCF files with bcftools (bcftools only)
    // [SPLIT from GLIMPSE_IMPUTE_CHUNKS]
    tag   "${chr}:${chunk_id}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(chr), val(chunk_id), path(bcf, arity: '1')

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
    tuple val(chr), val(chunk_id),
          path(bcf, arity: '1'),
          path("${bcf}.csi", arity: '1'),                                                                    emit: imputed
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process GLIMPSE_LIGATE {
    // Per-chr: ligate imputed chunks → one merged BCF per chromosome (GLIMPSE only)
    // Source: 01_mapping/3_LowPass_imputation/scripts/04_LigatePerChromsome.sh
    // [SPLIT] bcftools view/index moved to BCF_VCF process
    tag   "${chr}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/glimpse-bio:1.1.1--7517c3cd783f15ab'

    input:
    tuple val(chr), path(chunk_bcfs, arity: '1..*'), path(chunk_csis, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Sort chunks by filename for correct ligation order.
    ## Chunk IDs are zero-padded (00, 01, ...) so lexicographic sort is correct.
    ls DA.${chr}.imputed.*.bcf | sort > list.${chr}.txt

    GLIMPSE_ligate \\
        --input list.${chr}.txt \\
        --output DA.${chr}.merged.bcf

    # The GLIMPSE_ligate_static call above makes an index for the bcf.
    """

    stub:
    """
    touch DA.${chr}.merged.bcf
    """

    output:
    tuple val(chr),
          path("DA.${chr}.merged.bcf"),                                                                         emit: merged_bcf
    tuple val("${task.process}"), val('GLIMPSE'), eval("GLIMPSE_ligate --help | sed '/Version/!d; s/.* : //'"), emit: versions_glimpse,  topic: versions
}


process BCF_VCF {
    // Per-chr: convert merged BCF → bgzipped VCF + index (bcftools only)
    // [SPLIT from GLIMPSE_LIGATE]
    tag   "${chr}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(chr), path(bcfs, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools view ${bcfs} -Oz -o DA.${chr}.merged.vcf.gz
    bcftools index -f DA.${chr}.merged.vcf.gz
    """

    stub:
    """
    touch DA.${chr}.merged.vcf.gz
    touch DA.${chr}.merged.vcf.gz.csi
    """

    output:
    tuple val(chr),
          path("DA.${chr}.merged.vcf.gz"),
          path("DA.${chr}.merged.vcf.gz.csi"),                                                               emit: merged
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process GLIMPSE_QC_FILTER {
    // Per-chr: filter on imputation INFO score and annotate variant IDs (bcftools only)
    // Source: 01_mapping/3_LowPass_imputation/scripts/05_QCfilter_imputed_chr.sh
    // Note! you will need to remove intermediate bcf and vcf files after this step.
    // The final output is DA.${chr}.merged.qc_modi.vcf.gz
    // Note: bcftools-only process — uses highpass_snpfiltering container
    tag   "${chr}"
    label 'process_low'
    container 'community.wave.seqera.io/library/highpass_snpfiltering:0f456ef177340518'

    input:
    tuple val(chr), path(merged_vcf), path(merged_vcf_idx)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    bcftools filter -e 'INFO/INFO<0.8' ${merged_vcf} -Oz -o DA.${chr}.merged.qc.vcf.gz
    bcftools annotate --set-id '%CHROM:%POS' DA.${chr}.merged.qc.vcf.gz -Oz -o DA.${chr}.merged.qc_modi.vcf.gz
    bcftools index -f DA.${chr}.merged.qc_modi.vcf.gz
    """

    stub:
    """
    touch DA.${chr}.merged.qc_modi.vcf.gz
    touch DA.${chr}.merged.qc_modi.vcf.gz.csi
    """

    output:
    tuple val(chr),
          path("DA.${chr}.merged.qc_modi.vcf.gz"),
          path("DA.${chr}.merged.qc_modi.vcf.gz.csi"),                                                       emit: qc_vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '/bcftools/!d; s/.* //'"), emit: versions_bcftools, topic: versions
}


process VCF_TO_PLINK {
    // Per-chr: convert imputed VCF to PLINK binary format
    // Source: 01_mapping/3_LowPass_imputation/scripts/06_vcf_to_plink.sh
    //
    // first create plink files from vcf for each chr
    tag   "${chr}"
    label 'process_low'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    tuple val(chr), path(qc_vcf), path(qc_vcf_idx)

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #first create plink files from vcf for each chr
    plink --dog \\
        --vcf ${qc_vcf} \\
        --keep-allele-order \\
        --double-id \\
        --out DA_IMP_GENCOVE_${chr} \\
        --make-bed
    # --double-id causes both family and within-family IDs to be set to the sample ID
    """

    stub:
    """
    touch DA_IMP_GENCOVE_${chr}.bed
    touch DA_IMP_GENCOVE_${chr}.bim
    touch DA_IMP_GENCOVE_${chr}.fam
    touch DA_IMP_GENCOVE_${chr}.log
    touch DA_IMP_GENCOVE_${chr}.nosex
    """

    output:
    tuple val(chr),
          path("DA_IMP_GENCOVE_${chr}.bed"),
          path("DA_IMP_GENCOVE_${chr}.bim"),
          path("DA_IMP_GENCOVE_${chr}.fam"),                                                         emit: plink
    tuple val("${task.process}"), val('plink'), eval("plink --version | sed 's/PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions
}


process MERGE_CHR_PLINK {
    // Cohort-level: merge all per-chr PLINK files into final dataset
    // Source: 01_mapping/3_LowPass_imputation/scripts/07_Merge_chr_plink.sh
    tag   "cohort"
    label 'process_medium'
    container 'community.wave.seqera.io/library/plink:1.90b7.7--feb905bcf3811778'

    input:
    path(plink_files, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    ## Build merge list (chr2-38; chr1 is the --bfile base)
    ## make this list according to your sample names: mylist_DA_IMP_GENCOVE_chr1-38.txt
    ## DA_IMP_GENCOVE_chr1
    ## DA_IMP_GENCOVE_chr2
    ## DA_IMP_GENCOVE_chr3 etc.
    > merge_list.txt
    for i in \$(seq 2 38); do
        echo "DA_IMP_GENCOVE_chr\${i}" >> merge_list.txt
    done

    plink --allow-no-sex \\
        --dog \\
        --geno 0.05 \\
        --keep-allele-order \\
        --bfile DA_IMP_GENCOVE_chr1 \\
        --merge-list merge_list.txt \\
        --make-bed \\
        --out DA_IMP_GENCOVE_ALLCHR
    """

    stub:
    """
    touch DA_IMP_GENCOVE_ALLCHR.bed
    touch DA_IMP_GENCOVE_ALLCHR.bim
    touch DA_IMP_GENCOVE_ALLCHR.fam
    touch DA_IMP_GENCOVE_ALLCHR.log
    """

    output:
    path("DA_IMP_GENCOVE_ALLCHR.bed"),                                                               emit: bed
    path("DA_IMP_GENCOVE_ALLCHR.bim"),                                                               emit: bim
    path("DA_IMP_GENCOVE_ALLCHR.fam"),                                                               emit: fam
    path("DA_IMP_GENCOVE_ALLCHR.log"),                                                               emit: log
    tuple val("${task.process}"), val('plink'), eval("plink --version | sed 's/PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions
}
