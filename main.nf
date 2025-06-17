//**************Cut&Run pipeline*********

#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// === Workflow definition ===
workflow {
    // Channel of sample names
    samples_ch = Channel.fromList(params.SAMPLES)

    // Paired FASTQ file tuples
    fastqs_ch = samples_ch.map { sample ->
        def fq1 = file("${params.DATA}/${sample}*_*1.fq.gz")
        def fq2 = file("${params.DATA}/${sample}*_*2.fq.gz")
        tuple(sample, fq1, fq2)
    }

    // Pipe through each process
    trimmed_ch     = fastqs_ch        | cutadapt
    aligned_ch     = trimmed_ch       | bowtie2
    postaligned_ch = aligned_ch       | postAlign
    pebeds_ch      = postaligned_ch   | pebeds
    lengths_ch     = pebeds_ch        | lengthFiles
    gzbed_ch       = pebeds_ch        | zipBeds
    unique_ch      = gzbed_ch         | uniqBeds
}

// === Process: cutadapt ===
process cutadapt {
    tag "$sample"
    publishDir "${params.RESULTS}/$sample", mode: 'copy'
    module 'cutadapt/4.2'

    input:
      tuple val(sample), path(r1), path(r2)

    output:
      tuple val(sample),
            path("${sample}_R1_trim1.fastq.gz"),
            path("${sample}_R2_trim1.fastq.gz"),
            path("${sample}_R1_trim2.fastq.gz"),
            path("${sample}_R2_trim2.fastq.gz")

    script:
    """
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \\
             -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \\
             -o ${sample}_R1_trim1.fastq.gz \\
             -p ${sample}_R2_trim1.fastq.gz \\
             $r1 $r2
    
    //if read length <= 130, no trimming
    cutadapt --pair-filter=both \\
             -l ${params.len} \\
             -m ${params.m_len} \\
             -o ${sample}_R1_trim2.fastq.gz \\
             -p ${sample}_R2_trim2.fastq.gz \\
             ${sample}_R1_trim1.fastq.gz \\
             ${sample}_R2_trim1.fastq.gz
    """
}

// === Process: bowtie2 ===
process bowtie2 {
    tag "$sample"
    publishDir "${params.RESULTS}/$sample", mode: 'copy'
    module 'bowtie2/2.5.0'

    input:
      tuple val(sample), path(r1), path(r2)

    output:
      tuple val(sample), path("${sample}.sam")

    script:
    """
    bowtie2 -x ${params.GENOME} \\
            -1 $r1 \\
            -2 $r2 \\
            -S ${sample}.sam
    """
}

// === Process: postAlign ===
process postAlign {
    tag "$sample"
    publishDir "${params.RESULTS}/$sample", mode: 'copy'
    module 'samtools'
    module 'bedtools'

    input:
      tuple val(sample), path("${sample}.sam")

    output:
      tuple val(sample), path("${sample}.bam"), path("${sample}.bed")

    script:
    """
    samtools sort -n ${sample}.sam -o ${sample}.bam
    samtools view -bf 0x2 ${sample}.bam \
      | bedtools bamtobed -i stdin -bedpe > ${sample}.bed
    rm ${sample}.sam
    """
}

// === Process: pebeds ===
//# Create bed files for paired-end reads. Create pe.bed files using perl script that has chr no., chr START, chr END, fragment lengths 
process pebeds {
    tag "$sample"
    publishDir "${params.RESULTS}/$sample", mode: 'copy'
    module 'perl'

    input:
      tuple val(sample), path("${sample}.bed")

    output:
      tuple val(sample), path("${sample}_pe.bed")

    script:
    """
    perl ${params.SCRIPTS}/bedpe2bed.pl ${sample}.bed 2000 > ${sample}_pe.bed
    """
}

// === Process: lengthFiles ===
process lengthFiles {
    tag "$sample"
    publishDir "${params.RESULTS}/$sample", mode: 'copy'

    input:
      tuple val(sample), path("${sample}_pe.bed")

    output:
      tuple val(sample), path("${sample}_braw.len"), path("${sample}_bnorm.len")

    script:
    """
    awk '{print \$4}' ${sample}_pe.bed | sort -n | uniq -c | awk '{print \$2,\$1}' > ${sample}_braw.len
    len = `cat ${sample}_pe.bed | wc -l`
    awk -v len=\$len '{print \$1,\$2/len}' ${sample}_braw.len > ${sample}_bnorm.len
    """
}

// === Process: zipBeds ===
process zipBeds {
    tag "$sample"
    publishDir "${params.RESULTS}/$sample", mode: 'copy'

    input:
      tuple val(sample), path("${sample}_pe.bed")

    output:
      tuple val(sample), path("${sample}_pe.bed.gz")

    script:
    """
    gzip -c ${sample}_pe.bed > ${sample}_pe.bed.gz
    """
}

// === Process: uniqBeds ===
process uniqBeds {
    tag "$sample"
    publishDir "${params.RESULTS}/$sample", mode: 'copy'
    module 'perl'

    input:
      tuple val(sample), path("${sample}_pe.bed.gz")

    output:
      tuple val(sample), path("${sample}_uniq_pe.bed.gz")

    script:
    """
    perl ${params.SCRIPTS}/uniq_STDOUT.pl ${sample}_pe.bed.gz | gzip > ${sample}_uniq_pe.bed.gz
    """
}
