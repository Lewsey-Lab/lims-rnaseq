#!/usr/bin/env nextflow

params.reads_path = "${projectDir}/data/raw/reads/*1.fastq.gz"
// params.genome = "${projectDir}/data/raw/genome/*.@(fna|fa)"

Channel
    // closure specifies the first 9 characters as the accession id
    .fromFilePairs(params.reads_path, size: 1) {file -> file.name.substring(0,9)}
    .view()
    .set { reads_ch }

process fastqc {
    tag "FastQC on ${acc_id}"
    publishDir "${projectDir}/data/interim/fastqc", mode: "copy"

    input:
        tuple val(acc_id), path(read) from reads_ch

    output:
        path "${acc_id}" into fastqc_ch

    script:

        """
        mkdir -p ${acc_id}
        fastqc -o ${acc_id} ${read}
        """
}

fastqc_ch.view()