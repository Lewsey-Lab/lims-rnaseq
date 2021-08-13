#!/usr/bin/env nextflow

// in DSL2, processes are chained into workflows, and input and output blocks
// don't explicitly specify channels. The processes are written like functions
// in the workflow block, taking inputs as arguments, and other processes
// outputs as .out attributes
nextflow.enable.dsl = 2

params.reads = "${projectDir}/data/raw/reads/*.fastq.gz"
params.genome = "${projectDir}/data/raw/genome/*.{fna,fa}"

println """
        R N A S E Q - N F   P I P E L I N E
        ===================================
        reads: ${params.reads}
        genome: ${params.genome}
        """
        .stripIndent()

// Closure definition used whenever we want to distinguish between PE and SE
// entries by checking the length of their read list
def pe_criteria = branchCriteria {
    se: it[1].size() == 1
        return it + false
    pe: it[1].size() == 2
        return it + true
}

def create_fastq_channels(it) {
    def acc_id = it[0]
    def read_paths = it[1]
    def meta = [:]

    meta.id = acc_id
    meta.strandness = 'R'

    if (read_paths.size() == 1) {
        meta.single_end = true
    }

    else if (read_paths.size() == 2) {
        meta.single_end = false
    }

    def array = [meta, read_paths]
    return array
}

Channel
    // closure (like lambda functions in python) specifies the stream as tuple
    // of accession no. (single val) and filepath (list of one or two paths)
    .fromFilePairs(params.reads, size: -1) { file -> file.name.split('_')[0] }
    // branch splits stream into single read and paired read sub-streams, based
    // on whether the file list (index 1) has one or two files. A boolean is
    // appended to keep track of this information.
    // .branch(pe_criteria)
    .map {create_fastq_channels(it)}
    // Defines the channel name
    .set { raw_reads }

Channel
    .fromPath(params.genome)
    .set { genome }


process FASTQC {
    tag "FastQC on ${acc_id}"
    // Copies outputs of process to publishdir. No need to use absolute paths
    // in shell/script block
    publishDir "${projectDir}/reports/fastqc", mode: 'copy'

    input:
        // Breaks the tuple stream to acc_id and read path
        tuple val(acc_id), path(reads), val(pe)

    output:
        // Folder containing fastqc html and zip. '*' is necessary as only the
        // contents of the folder is unique, otherwise fastqc folders of the
        // same name from the trimmed output would result in input collisions
        path "${acc_id}/*"

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

    // --fullnames arg used otherwise multiqc will not distinguish untrimmed
    // and trimmed reads
    shell:
        '''
        multiqc --fullnames .
        '''
}

process TRIMMING {
    tag "${acc_id}"

    // Separates the trimmed reads from the generated reports
    publishDir(
        path: "${projectDir}/data/interim/trim_reads",
        pattern: '*_trimmed.fq.gz',
        mode: 'copy'
        )

    publishDir(
        path:"${projectDir}/reports/trim_galore",
        pattern: '*.txt',
        mode: 'copy'
    )

    input:
        tuple val(acc_id), path(reads), val(pe)

    output:
        tuple val(acc_id), path('*.fq.gz'), val(pe), emit: trim_reads
        path('*.txt')

    // TODO: Provide core arguments dynamically
    shell:
        if ( pe == true ) {
            '''
            trim_galore --illumina --trim-n --cores 2 --paired !{reads[0]} !{reads[1]}
            '''
        }

        else {
            '''
            trim_galore --illumina --trim-n --cores 2 !{reads[0]}
            '''
        }
}

process INDEXING {
    tag "${genome.name}"

    storeDir "${projectDir}/data/interim/index/${genome.name}"

    input:
        path(genome)

    output:
        path('*')

    // TODO: provide core argument dynamically
    shell:
        '''
        hisat2-build -p 12 !{genome} !{genome}
        '''
}

process ALIGNMENT {
    tag "${acc_id}"

    publishDir(
        path: "${projectDir}/data/processed/hisat2",
        pattern: "*.ba[m,i]",
        mode: 'copy'
    )

    publishDir(
        path: "${projectDir}/reports/hisat2",
        pattern: "*.txt",
        mode: 'copy'
    )

    input:
        tuple val(acc_id), path(reads), val(pe)
        path(index)

    output:
        path('*.txt')
        path('*.bam')
        path('*bai')

    shell:
        index_name = index[0].name.split(/.\d+.ht2/)[0]
        if (pe == true) {
            arg = "-1 ${reads[0]} -2 ${reads[1]}"
        }
        else {
            arg = "-U ${reads[0]}"
        }
        '''
        hisat2 \
        --min-intronlen 20 --max-intronlen 6000 \
        --rna-strandness R \
        -p 10 \
        --new-summary --summary-file !{acc_id}.txt \
        -x !{index_name} \
        !{arg} \
        | samtools sort -O BAM \
        |tee !{acc_id}.sorted.max.intron.6000.bam \
        | samtools index - !{acc_id}.sorted.max.intron.6000.bai
        '''
}

workflow {
    // raw_reads.se and raw_reads.pe are mixed as they're processed the same
    raw_reads.view()
    // TRIMMING(raw_reads.mix())
    // // Mixing trimmed and untrimmed reads before fastqc
    // all_reads = TRIMMING.out.trim_reads.mix(raw_reads.se, raw_reads.pe)
    // FASTQC(all_reads)
    // MULTIQC(FASTQC.out.collect())
    // INDEXING(genome)
    // ALIGNMENT(TRIMMING.out.trim_reads, INDEXING.out.first())
}
