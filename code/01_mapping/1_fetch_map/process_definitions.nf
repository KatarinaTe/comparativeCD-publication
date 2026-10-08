// ============================================================================
// process_definitions.nf — 01_mapping / 1_fetch_map
//
// Per-run/per-sample processes: FASTQ -> sorted BAM -> merged/BQSR'd BAM + gVCF
// -> per-sample stats. Runs identically for HighPass and LowPass samples — track
// distinction is which samplesheet/invocation calls this pipeline, not separate
// process aliases (see 1_fetch_map.nf).
//
// Split out of the original code/01_mapping/process_definitions.nf. Process script
// bodies are unchanged from that file — only which pipeline file they live in changed.
// See archive/conversion_notes/01_mapping.md for the split rationale.
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


process DOG10K_MAPPING {
    // Per-run: FASTQ → sorted BAM
    // Source: 01_mapping/1_High-LowPass_mapping/dog10k_mapping.sh (mapping section)
    tag   "${sample_id}:${run_id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89'

    input:
    tuple val(sample_id), val(run_id), path(r1_path, arity: '1'), path(r2_path, arity: '1')
    tuple path(assembly_ref), path(assembly_ref_idx, arity: '1..*')

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${sample_id}_${run_id}"
    """
    ## mapping
    bwa-mem2 mem \\
        -K 100000000 \\
        -t ${task.cpus} \\
        -Y \\
        -R '@RG\\tID:${run_id}\\tPL:illumina\\tLB:${sample_id}\\tSM:${sample_id}' \\
        ${assembly_ref} \\
        ${r1_path} \\
        ${r2_path} \\
        | samtools view -hbS > ${prefix}.bam

    samtools sort \\
        -@ ${task.cpus} \\
        -m 1G \\
        ${prefix}.bam \\
        -o ${prefix}.sorted.bam
    rm -f ${prefix}.bam
    samtools index ${prefix}.sorted.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${sample_id}_${run_id}"
    """
    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    """

    output:
    tuple val(sample_id),
          path("${sample_id}_${run_id}.sorted.bam",     arity: '1'),
          path("${sample_id}_${run_id}.sorted.bam.bai", arity: '1'),                                  emit: bam
    tuple val("${task.process}"), val('bwa-mem2'),  eval("bwa-mem2 version"),                         emit: versions_bwamem2,  topic: versions
    tuple val("${task.process}"), val('samtools'),  eval("samtools --version | sed '1!d; s/.* //'"),  emit: versions_samtools, topic: versions
}


process DOG10K_MERGE_MARKDUPS_BQSR_GVCF {
    // Per-sample: merge all runs → MarkDups → BQSR → HaplotypeCaller → gVCF
    // Source: 01_mapping/1_High-LowPass_mapping/dog10k_mapping.sh (merge/BQSR/GVCF)
    tag   "${sample_id}"
    label 'process_high'
    container 'community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89'

    input:
    tuple val(sample_id), path(bams, arity: '1..*'), path(bais, arity: '1..*')
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
    ## [FIX 2026-09-19] The original script (and this port, faithfully, until now) referenced a
    ## chunk_11.g.vcf.gz in the GatherVcfs call below that nothing ever produced -- chunks_dir only
    ## ever has 10 bed files, and data/raw-data/references/chunks_10/10chunk (the manifest those 10
    ## files were generated from) assigns all 43 real contigs (chr1-38, chrX, both chrY variants,
    ## chrM) to buckets 1-10 only; no bucket 11 exists in it. The thousands of small chrUn_* unplaced
    ## scaffolds in the real assembly .fai were never assigned to any bucket either, chunk_11
    ## included -- they were out of scope for the whole 10-chunk design from the start, not
    ## specifically missing from an unfinished 11th one. Removed the dangling
    ## `-I ...chunk_11.g.vcf.gz` line below; nothing is lost, since every one of the 43 real contigs
    ## is already gathered via the 43 explicit -I ...sc.<contig>.g.vcf.gz lines. See
    ## REPRODUCIBILITY_AUDIT.md's 2026-09-19 update for the full investigation.
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
    tuple val(sample_id),
          path("*.sorted.merged.MarkDups.BQSR.bam"),
          path("*.sorted.merged.MarkDups.BQSR.bam.bai"),                                                                      emit: bam
    tuple val(sample_id),
          path("*.sorted.merged.MarkDups.BQSR.g.vcf.gz"),
          path("*.sorted.merged.MarkDups.BQSR.g.vcf.gz.tbi"),                                                                 emit: gvcf
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/.* //'"),                           emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('picard'),   eval("picard MarkDuplicates --version 2>&1 | sed '/Version/!d; s/.*://'"), emit: versions_picard,   topic: versions
    tuple val("${task.process}"), val('gatk'),     eval("gatk --version | sed '/GATK/!d;s/.* v//'"),                          emit: versions_gatk,     topic: versions
    tuple val("${task.process}"), val('tabix'),    eval("tabix --version | sed '/tabix/!d; s/.* //'"),                        emit: versions_tabix,    topic: versions
}


process RUN_STATS {
    // Per-sample: BAM/gVCF quality metrics and depth summary
    // Source: 01_mapping/1_High-LowPass_mapping/run-stats.py
    //
    // summarize_stats() is not functional:
    //   it reads MarkDups metrics files whose source is unknown (not documented
    //   in the original scripts). The call is commented out below; the rest of
    //   summarize_stats() (file sizes, depth stats, X/Auto ratio) runs correctly.
    tag   "${sample_id}"
    label 'process_medium'
    container 'community.wave.seqera.io/library/dog10k_mapping:f7e5c65de8fdfe89'

    input:
    tuple val(sample_id), path(bam, arity: '1'), path(bai, arity: '1'),
                          path(gvcf, arity: '1'), path(gvcf_tbi, arity: '1')
    tuple path(assembly_ref), path(assembly_ref_idx, arity: '1..*')
    tuple path(depth_sites),  path(depth_sites_tbi, arity: '1')

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 - << 'END_PYTHON_SCRIPT'
    # program to process stats from CRAM file
    # Source: 01_mapping/1_High-LowPass_mapping/run-stats.py

    import sys
    import subprocess
    import os
    import argparse
    import time
    import socket
    import shutil
    import gzip
    import numpy as np

    ###############################################################################
    # Helper function to run commands, handle return values and print to log file
    def runCMD(cmd):
        val = subprocess.Popen(cmd, shell=True).wait()
        if val == 0:
            pass
        else:
            print('command failed')
            print(cmd)
            sys.exit(1)
    ###############################################################################
    # setup paths to default programs to use and checks for required programs
    def check_prog_paths(myData):
        print('checking required programs')
        for p in ['gatk','samtools','Rscript','tabix']:
            if shutil.which(p) is None:
                s = p + ' not found in path! please fix (module load?)'
                print(s, flush=True)
                sys.exit()
            else:
                print('%s\\t%s' % (p,shutil.which(p)),flush=True)
    #############################################################################
    def run_flagstat(myData):
        print('\\nrunning flagstat ...', flush=True)
        myData['flagStatFileName'] = myData['cramFileName'] + '.flagstat'

        if os.path.isfile(myData['flagStatFileName']) is False:
            cmd = 'samtools flagstat %s > %s ' % (myData['cramFileName'],myData['flagStatFileName'])
            print(cmd,flush=True)
            runCMD(cmd)
            print('Done',flush=True)
        else:
            print('%s exists, skipping flagstat' % myData['flagStatFileName'],flush=True)
    ###############################################################################
    def run_insertmetrics(myData):
        print('\\nrunning CollectInsertSizeMetrics ...', flush=True)

        myData['sizeMetricsFileName'] = myData['cramFileName'] + '.insert_size_metrics.txt'
        myData['sizeMetricsFileHistName'] = myData['cramFileName'] + '.insert_size_metrics.hist.pdf'

        if os.path.isfile(myData['sizeMetricsFileName']) is False:
            cmd = 'gatk CollectInsertSizeMetrics -R %s -I %s -O %s -H %s' % (myData['ref'],myData['cramFileName'],
                                                                             myData['sizeMetricsFileName'],myData['sizeMetricsFileHistName'] )
            print(cmd,flush=True)
            runCMD(cmd)
            print('Done',flush=True)
        else:
            print('%s exists, skipping CollectInsertSizeMetrics' % myData['sizeMetricsFileName'],flush=True)
    ###############################################################################
    def run_alignmentmetrics(myData):
        print('\\nrunning CollectAlignmentSummaryMetrics ...', flush=True)

        myData['alignMetricsFileName'] = myData['cramFileName'] + '.alignment_summary_metrics'

        if os.path.isfile(myData['alignMetricsFileName']) is False:
            cmd = 'gatk CollectAlignmentSummaryMetrics -R %s -I %s -O %s ' % (myData['ref'],myData['cramFileName'],
                                                                             myData['alignMetricsFileName'], )
            print(cmd,flush=True)
            runCMD(cmd)
            print('Done',flush=True)
        else:
            print('%s exists, skipping CollectAlignmentSummaryMetrics' % myData['alignMetricsFileName'],flush=True)
    ###############################################################################
    def run_index_cram(myData):
        print('\\nrunning index cram ...', flush=True)
        craiFileName = myData['cramFileName'] + '.crai'
        if os.path.isfile(craiFileName) is False:
            cmd = 'samtools index %s' % myData['cramFileName']
            print(cmd,flush=True)
            runCMD(cmd)
            print('Done',flush=True)
        else:
            print('%s exists, skipping index cram' % craiFileName,flush=True)
    ###############################################################################
    def run_index_gvcf(myData):
        print('\\nrunning index gvcf ...', flush=True)
        tbiFileName = myData['gvcf'] + '.tbi'

        if os.path.isfile(tbiFileName) is False:
            cmd = 'tabix -p vcf %s' %  myData['gvcf']
            print(cmd,flush=True)
            runCMD(cmd)
            print('Done',flush=True)
        else:
            print('%s exists, skipping index gvcf' % tbiFileName,flush=True)
    ###############################################################################
    def run_genotype_known_sites(myData):
        print('\\nrunning genotype known sites ...', flush=True)
        myData['knownSitesVCF'] = myData['cramFileName'] + '.knownsites.vcf.gz'

        if os.path.isfile(myData['knownSitesVCF']) is False:
            cmd = 'gatk --java-options "-Xmx2g" GenotypeGVCFs '
            cmd += '-R %s -V %s -O %s ' % (myData['ref'],myData['gvcf'],myData['knownSitesVCF'] )
            cmd += ' --include-non-variant-sites --intervals %s ' % (myData['sites'])
            print(cmd,flush=True)
            runCMD(cmd)
            print('Done',flush=True)
        else:
            print('%s exists, genotype known sites' % myData['knownSitesVCF'],flush=True)
    ###############################################################################
    def calc_effective_depth(myData):
        print('\\nrunning calc effective depth ...', flush=True)
        myData['depthSummary'] = myData['cramFileName'] + '.knownsites.depth.txt'

        if os.path.isfile(myData['depthSummary']) is True:
            print('%s exists' % myData['depthSummary'],flush=True)
            return
        autoDp = []
        xDp = []
        inFile = gzip.open(myData['knownSitesVCF'],'rt')
        for line in inFile:
            if line[0] == '#':
                continue
            line = line.rstrip()
            line = line.split()

            infoField = line[7]
            infoField = infoField.split(';')
            dp = -1
            for i in infoField:
                if i[0:3] == 'DP=':
                    dp = int(i.split('=')[1])
                    break

            if dp == -1:
                dp = 0
            if line[0] == 'chrX':
                xDp.append(dp)
            else:
                autoDp.append(dp)
        inFile.close()

        outFile = open(myData['depthSummary'],'w')

        outFile.write('#chrom\\ttotalSites\\tMean\\tStd\\tMedian\\n')
        outFile.write('Autos\\t%i\\t%.2f\\t%.2f\\t%.1f\\n' % (len(autoDp),np.mean(autoDp),np.std(autoDp),np.median(autoDp) ))
        outFile.write('ChrX\\t%i\\t%.2f\\t%.2f\\t%.1f\\n' % (len(xDp),np.mean(xDp),np.std(xDp),np.median(xDp) ))
        outFile.close()
    ###############################################################################
    def summarize_stats(myData):
        myData['statsSummary'] = myData['cramFileName'] + '.stats.txt'

        outFile = open(myData['statsSummary'],'w')

        sn = myData['cramFileName'].split('/')[-1].split('.')[0]

        outFile.write('SampleName\\t%s\\n' % (sn))

        cramFileSize = os.path.getsize(myData['cramFileName']) / (1024*1024*1024)

        outFile.write('CramSize\\t%.2f Gb\\n' % cramFileSize)

        cramFileSize = os.path.getsize(myData['gvcf']) / (1024*1024*1024)
        outFile.write('GVCFSize\\t%.2f Gb\\n' % cramFileSize)

        # now get coverage stats
        inFile = open(myData['depthSummary'],'r')
        lines=inFile.readlines()
        inFile.close()

        autoLine = lines[1].rstrip().split()
        xLine = lines[2].rstrip().split()

        outFile.write('effectiveAutoMean\\t%s\\neffectiveAutoMedian\\t%s\\n' % (autoLine[2],autoLine[4]) )
        outFile.write('effectiveXMean\\t%s\\neffectiveXMedian\\t%s\\n' % (xLine[2],xLine[4]) )

        xvsAutoMean =  float(xLine[2]) / float(autoLine[2])

        if float(autoLine[4]) == 0.0:
            xvsAutoMedian = 0.0
        else:
            xvsAutoMedian = float(xLine[4]) / float(autoLine[4])

        outFile.write('Mean(X/Auto)\\t%.2f\\n' % xvsAutoMean)
        outFile.write('Median(X/Auto)\\t%.2f\\n' % xvsAutoMedian)

        # ── Non-functional: reads MarkDups metrics files whose source is unknown
        #    (not documented in the original scripts). The entire section after
        #    the X/Auto ratio is non-functional because it depends on these files.

        # dupMetFileName = myData['cramFileName'].replace('.cram','.sort.md.metricts.txt')
        # inFile = open(dupMetFileName,'r')
        # lines = inFile.readlines()
        # inFile.close()
        # metLine = lines[7].rstrip().split()
        # dupF = metLine[8]
        # outFile.write('fractionDup\\t%s\\n' % dupF)

        # inFile = open(myData['alignMetricsFileName'],'r')
        # lines = inFile.readlines()
        # pairLine = lines[9].rstrip().split()
        # inFile.close()
        # outFile.write('totalPairedReads\\t%s\\n' % pairLine[1])
        # outFile.write('fractionAligned\\t%s\\n' % pairLine[6])
        # outFile.write('mismatchRate\\t%s\\n' % pairLine[12])
        # outFile.write('indelRate\\t%s\\n' % pairLine[14])
        # outFile.write('meanReadLen\\t%s\\n' % pairLine[15])
        # outFile.write('fractionImproperPairs\\t%s\\n' % pairLine[19])
        # outFile.write('fractionChimera\\t%s\\n' % pairLine[22])

        # inFile = open(myData['sizeMetricsFileName'],'r')
        # lines = inFile.readlines()
        # inFile.close()
        # statLine = lines[7].rstrip().split()
        # numPairsFirstLine= int(statLine[7])
        # totPairs = numPairsFirstLine
        # i = 8
        # while lines[i] != '\\n':
        #    nl = lines[i].rstrip().split()
        #    totPairs += int(nl[7])
        #    i += 1
        # fractionPairsAssigned = numPairsFirstLine/totPairs
        # outFile.write('pairOrientation\\t%s\\n' % statLine[8])
        # outFile.write('fractionWithPairOrientation\\t%.4f\\n' % fractionPairsAssigned )
        # outFile.write('meanInsertLen\\t%s\\n' % statLine[5])
        # outFile.write('stdInsertLen\\t%s\\n' % statLine[6])
        # outFile.write('medianInsertLen\\t%s\\n' % statLine[0])
        # outFile.write('madInsertLen\\t%s\\n' % statLine[2])
        # outFile.close()

        print('Summary written to',myData['statsSummary'])
    ###############################################################################

    # SETUP — paths are injected by Nextflow interpolation (replaces argparse)
    myData = {} # dictionary for keeping and passing information
    myData['cramFileName'] = '${bam}'
    myData['ref'] = '${assembly_ref}'
    myData['gvcf'] = '${gvcf}'
    myData['sites'] = '${depth_sites}'

    # make sure programs are available
    check_prog_paths(myData)
    run_flagstat(myData)
    run_insertmetrics(myData)
    run_alignmentmetrics(myData)
    run_index_cram(myData)
    run_index_gvcf(myData)
    run_genotype_known_sites(myData)
    calc_effective_depth(myData)
    # summarize_stats(myData)  # commented out — non-functional (see process comment)
    END_PYTHON_SCRIPT
    """

    stub:
    """
    touch ${bam}.flagstat
    touch ${bam}.insert_size_metrics.txt
    touch ${bam}.insert_size_metrics.hist.pdf
    touch ${bam}.alignment_summary_metrics
    touch ${bam}.knownsites.vcf.gz
    touch ${bam}.knownsites.vcf.gz.tbi
    touch ${bam}.knownsites.depth.txt
    # touch ${bam}.stats.txt  // commented out — summarize_stats() produces this file (non-functional, see process comment)
    """

    output:
    tuple val(sample_id), path("*.bam.knownsites.depth.txt",         arity: '1'),                            emit: depth
    // tuple val(sample_id), path("*.bam.stats.txt",                    arity: '1'),                            emit: stats  // commented out — summarize_stats() is non-functional
    tuple val(sample_id), path("*.bam.flagstat",                     arity: '1'),                            emit: flagstat
    tuple val(sample_id), path("*.bam.insert_size_metrics.txt",      arity: '1'),                            emit: insert_metrics
    tuple val(sample_id), path("*.bam.insert_size_metrics.hist.pdf", arity: '1'),                            emit: insert_hist
    tuple val(sample_id), path("*.bam.alignment_summary_metrics",    arity: '1'),                            emit: align_metrics
    tuple val(sample_id), path("*.bam.knownsites.vcf.gz"),
                          path("*.bam.knownsites.vcf.gz.tbi"),                                               emit: knownsites_vcf
    tuple val("${task.process}"), val('python3'),   eval("python3 --version | sed '/Python/!d; s/.* //'"),   emit: versions_python,   topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed '1!d; s/.* //'"),          emit: versions_samtools, topic: versions
    tuple val("${task.process}"), val('gatk'),     eval("gatk --version | sed '/GATK/!d;s/.* v//'"),         emit: versions_gatk,     topic: versions
    tuple val("${task.process}"), val('tabix'),    eval("tabix --version | sed '/tabix/!d; s/.* //'"),       emit: versions_tabix,    topic: versions
    tuple val("${task.process}"), val('rscript'),  eval("Rscript --version | sed 's/.*version //;s/ .*//'"), emit: versions_rscript,  topic: versions    
}


