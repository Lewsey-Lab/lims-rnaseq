#!/usr/bin/env nextflow

// in DSL2, processes are chained into workflows, and input and output blocks
// don't explicitly specify channels. The processes are written like functions
// in the workflow block, taking inputs as arguments, and other processes
// outputs as .out attributes
nextflow.enable.dsl = 2

params.reads = "${projectDir}/data/raw/reads/*.fastq.gz"
// params.genome = "${projectDir}/data/raw/genome/*.@(fna|fa)"

println """
        R N A S E Q - N F   P I P E L I N E
        ===================================
        reads: ${params.reads}
        """
        .stripIndent()

Channel
    // closure (like lambda functions in python) specifies the stream as tuple
    // of accession no. (single val) and filepath (list of one or two paths)
    .fromFilePairs(params.reads, size: -1) { file -> file.name.substring(0, 9) }
    // branch splits stream into single read and paired read sub-streams, based
    // on whether the file list (index 1) has one or two files. A boolean is
    // appended to keep track of this information.
    .branch {
        se: it[1].size() == 1
            return it + false
        pe: it[1].size() == 2
            return it + true
    }
    // Defines the channel name
    .set { reads_ch }

process FASTQC {
    tag "FastQC on ${acc_id}"
    // Copies outputs of process to publishdir. No need to use absolute paths
    // in shell/script block
    publishDir "${projectDir}/reports/fastqc", mode: 'copy'

    input:
        // Breaks the tuple stream to acc_id and read path
        tuple val(acc_id), path(reads)

    output:
        // Folder containing fastqc html and zip
        path "${acc_id}"

    // shell/script block does work in 'work' temp directory, NOT the project
    // directory (important to note!) hence the folder passed to fastqc as
    // output dir only needs to be one level deep. This is then 'published' by
    // publishDir, no need to muck around with absolute/relative directories in
    // the shell/script block
    shell:
        '''
        mkdir -p !{acc_id}
        fastqc -q -t 6 -o !{acc_id} !{reads}
        '''
}

process MULTIQC {
    publishDir "${projectDir}/reports/fastqc", mode: 'copy'

    input:
        path('*')

    output:
        path('multiqc_report.html')

    shell:
        '''
        multiqc .
        '''
}

process TRIMMING {
    tag "${acc_id}"

    publishDir(
        path: "${projectDir}/data/interim/trim_reads",
        pattern: '*_trimmed.fastq.gz',
        mode: 'copy'
        )

    publishDir(
        path:"${projectDir}/reports/trim_galore",
        pattern: '*.txt',
        mode: 'copy'
    )

    input:
        tuple val(acc_id), path(reads), val(pe)

    // output:
    //     path('*_trimmed.fastq.gz'), emit: trim_reads
    //     path('*.txt')

    shell:
        println "acc_id: ${acc_id} path: ${reads}, PE ${pe}"
        '''echo test'''
// if ( pe == true )
//     '''
//     trim_galore --illumina --trim-n !{reads[0]} !{reads[1]}
//     '''
// else
//     '''
//     trim_galore --illumina --trim-n !{reads[0]}
//     '''
}

workflow {
    reads_ch.mix().view()
    // reads_ch.se and reads_ch.pe are mixed as they're processed the same
    // FASTQC(reads_ch.mix())
    // MULTIQC(FASTQC.out.collect())
    TRIMMING(reads_ch.mix())
}
