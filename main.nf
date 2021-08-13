#!/usr/bin/env nextflow

// in DSL2, processes are chained into workflows, and input and output blocks
// don't explicitly specify channels. The processes are written like functions
// in the workflow block, taking inputs as arguments, and other processes
// outputs as .out attributes
nextflow.enable.dsl = 2

params.reads = "${projectDir}/data/raw/reads/*.fastq.gz"
params.genome = "${projectDir}/data/raw/genome/*.{fna,fa}"
params.dev = false
params.number_of_inputs = 2

println """
        R N A S E Q - N F   P I P E L I N E
        ===================================
        reads: ${params.reads}
        genome: ${params.genome}
        """
        .stripIndent()

def create_metadata(it) {
    // Infers single-end reads and combines it with ID in a meta dict
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
    // uses preexisting .fromFilePairs method to get [acc_id, [reads]]
    .fromFilePairs(params.reads, size: -1) {file -> file.name.split('_')[0]}
    // if params.dev is false, take all inputs
    .take(params.dev ? params.number_of_inputs : -1)
    // create metadata and return [meta, [reads]]
    .map{create_metadata(it)}
    .set{raw_reads}

Channel
    .fromPath(params.genome)
    // Takes only the first genome it finds in param.genome
    .first()
    .set{genome}


process FASTQC {
    tag "${meta.id}"
    // Copies outputs of process to publishdir. No need to use absolute paths
    // in shell/script block
    publishDir "${projectDir}/reports/fastqc", mode: 'copy'

    input:
        // Breaks the tuple stream to meta dict and read path
        tuple val(meta), path(reads)

    output:
        // Folder containing fastqc html and zip. '*' is necessary as only the
        // contents of the folder is unique, otherwise fastqc folders of the
        // same name from the trimmed output would result in input collisions
        path "${meta.id}/*"

    // shell/script block does work in 'work' temp directory, NOT the project
    // directory (important to note!) hence the folder passed to fastqc as
    // output dir only needs to be one level deep. This is then 'published' by
    // publishDir, no need to muck around with absolute/relative directories in
    // the shell/script block
    shell:
        '''
        mkdir -p !{meta.id}
        fastqc -q -t 6 -o !{meta.id} !{reads}
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
    tag "${meta.id}"

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
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path('*.fq.gz'), emit: trim_reads
        path('*.txt')

    // TODO: Provide core arguments dynamically
    shell:
        if ( meta.single_end == true ) {
            args = "${reads[0]}"
        }
        else {
            args = "--paired ${reads[0]} ${reads[1]}"
        }
        '''
        trim_galore --illumina --trim-n --cores 2 !{args}
        '''
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
    tag "${meta.id}"

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
        tuple val(meta), path(reads)
        path(index)

    output:
        path('*.txt')
        path('*.bam')
        path('*bai')

    shell:
        index_name = index[0].name.split(/.\d+.ht2/)[0]
        if (meta.single_end == true) {
            arg = "-U ${reads[0]}"
        }
        else {
            arg = "-1 ${reads[0]} -2 ${reads[1]}"
        }
        '''
        hisat2 \
        --min-intronlen 20 --max-intronlen 6000 \
        --rna-strandness R \
        -p 10 \
        --new-summary --summary-file !{meta.id}.txt \
        -x !{index_name} \
        !{arg} \
        | samtools sort -O BAM \
        |tee !{meta.id}.sorted.max.intron.6000.bam \
        | samtools index - !{meta.id}.sorted.max.intron.6000.bai
        '''
}

workflow {
    TRIMMING(raw_reads)
    // // Mixing trimmed and untrimmed reads before fastqc
    all_reads = raw_reads.mix(TRIMMING.out.trim_reads)
    FASTQC(all_reads)
    MULTIQC(FASTQC.out.collect())
    INDEXING(genome)
    ALIGNMENT(TRIMMING.out.trim_reads, INDEXING.out)
}
