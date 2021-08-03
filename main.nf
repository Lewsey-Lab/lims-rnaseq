#!/usr/bin/env nextflow

// in DSL2, processes are chained into workflows, and input and output blocks
// don't explicitly specify channels. The processes are written like functions
// in the workflow block, taking inputs as arguments, and other processes
// outputs as .out attributes
nextflow.enable.dsl=2

params.reads = "${projectDir}/data/raw/reads/*.fastq.gz"
// params.genome = "${projectDir}/data/raw/genome/*.@(fna|fa)"

println """
        R N A S E Q - N F   P I P E L I N E
        ===================================
        reads: ${params.reads}
        """
        .stripIndent()

Channel
    // closure specifies the stream as tuple of accession no. and filepath
    .fromFilePairs(params.reads, size: 1) {file -> file.name.substring(0,9)}
    // Defines the channel name
    .set { reads_ch }

process FASTQC {
    tag "FastQC on ${acc_id}"
    // Copies outputs of process to publishdir. No need to use absolute paths
    // in shell/script block
    publishDir "${projectDir}/reports/fastqc", mode: "copy"

    input:
        // Breaks the tuple stream to acc_id and read path
        tuple val(acc_id), path(read)

    output:
        // Folder containing fastqc html and zip
        path "${acc_id}"

    // shell/script block does work in 'work' temp directory, NOT the project
    // directory (important to note!) hence the folder passed to fastqc as
    // output dir only needs to be one level deep. This is then 'published' by
    // publishDir, no need to muck around with absolute/relative directories in
    // the shell/script block
    shell:
        """
        mkdir -p !{acc_id}
        fastqc -q -t 6 -o !{acc_id} !{read}
        """
}

process MULTIQC {
    publishDir "${projectDir}/reports/fastqc", mode: "copy"

    input:
        path "*"

    output:
        path "multiqc_report.html"

    shell:
    """
    multiqc .
    """
}

workflow {
    FASTQC(reads_ch)
    MULTIQC(FASTQC.out.collect())
}