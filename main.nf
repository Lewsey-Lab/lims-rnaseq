#!/usr/bin/env nextflow

params.reads_path = "${projectDir}/data/raw/reads/*.fastq.gz"
// params.genome = "${projectDir}/data/raw/genome/*.@(fna|fa)"

Channel
    // closure specifies the first 9 characters as the accession id
    .fromFilePairs(params.reads_path, size: 1) {file -> file.name.substring(0,9)}
    .view()
    .set { reads_ch }

process fastqc {
    tag "FastQC on ${acc_id}"
    // publishDir "${projectDir}/data/interim/fastqc"

    input:
        tuple val(acc_id), path(read) from reads_ch

    // output:
    //     path "${projectDir}/data/interim/fastqc/${acc_id}" into fastqc_ch

    script:

        outdir = "${projectDir}/data/interim/fastqc/${acc_id}"

        """
        mkdir -p ${outdir}
        fastqc -o ${outdir} ${read}
        """
}
