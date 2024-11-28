#!/usr/bin/env nextflow

params.input = 'samplesheet.csv'
params.fasta = '/camp/lab/bauerd/home/users/chakraa2/projects/harriet/202405/ref/amplicon.fa'

// Create channel for input samplesheet
Channel
    .fromPath(params.input)
    .splitCsv(header:true)
    .map { row -> [ row.sample, file(row.fastq1, checkIfExists: true), file(row.fastq2, checkIfExists: true) ] }
    .set { ch_input }

Channel
    .fromPath(params.fasta)
    .set { ch_fasta }

process FASTP { 

    tag "${sample}"
    label "process_medium"
    container 'quay.io/biocontainers/fastp:0.23.4--h5f740d0_0'

    input: 
    tuple val(sample), path(fastq1), path(fastq2)

    output: 
    tuple val(sample), path("${sample}.trimmed.R1.fq.gz"), path("${sample}.trimmed.R2.fq.gz")

    script: 
    """
    fastp \
    -i ${fastq1} \
    -I ${fastq2} \
    -o ${sample}.trimmed.R1.fq.gz \
    -O ${sample}.trimmed.R2.fq.gz
    """
} 

process BOWTIE2_INDEX {

    tag "$fasta"
    label "process_medium"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    input:
    path(fasta)

    output:
    path("${fasta.simpleName}.*.bt2")

    script:
    """
    bowtie2-build $fasta ${fasta.simpleName}
    """
}

process BOWTIE2_ALIGN { 

    tag "${sample}"
    label "process_medium"
    publishDir "${params.outdir}/bowtie2", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    input: 
    tuple val(sample), path(fastq1), path(fastq2)
    path(index)

    output: 
    tuple val(sample), path("${sample}.bam"), path("${sample}.bam.bai")

    script: 
    """
    bowtie2 \
    -p ${task.cpus} \
    -k 1 --no-unal \
    -x ${index[0].simpleName} \
    -1 ${fastq1}  \
    -2 ${fastq2} | \
    samtools view -hb /dev/stdin | \
    samtools sort /dev/stdin > ${sample}.bam && \
    samtools index ${sample}.bam
    """
} 

process SAMTOOLS_IDXSTATS {

    tag "${sample}"
    label "process_medium"
    publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    input: 
    tuple val(sample), path(bam), path(bai)

    output: 
    tuple val(sample), path("${sample}.idxstats.txt")

    script: 
    """
    samtools idxstats ${bam} > ${sample}.idxstats
    awk '{{OFS="\t"}}{{print "${sample}", \$1, \$3}}' ${sample}.idxstats > ${sample}.idxstats.txt
    """
}

process SUMMARISE_IDXSTATS {

    tag "idxstats"
    label "process_low"
    publishDir "${params.outdir}/samtools", mode: 'copy', overwrite: true
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    input: 
    path(idxstats)

    output: 
    path("summary.idxstats.txt")

    script: 
    """
    cat $idxstats | sort -k 1,1 > summary.idxstats.txt
    """
}

workflow {

    FASTP(ch_input)
    BOWTIE2_INDEX(ch_fasta)
    BOWTIE2_ALIGN(FASTP.out, BOWTIE2_INDEX.out.collect())
    SAMTOOLS_IDXSTATS(BOWTIE2_ALIGN.out)
    SUMMARISE_IDXSTATS(SAMTOOLS_IDXSTATS.out.map{it[1]}.collect())

}