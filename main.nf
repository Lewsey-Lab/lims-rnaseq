#!/usr/bin/env nextflow

// in DSL2, processes are chained into workflows, and input and output blocks
// don't explicitly specify channels. The processes are written like functions
// in the workflow block, taking inputs as arguments, and other processes
// outputs as .out attributes
nextflow.enable.dsl = 2

params.genome = false
if (!params.genome) {exit 1, "Please specify genome folder name with --genome <genome_name>"}

params.genomeDir = "${launchDir}/data/genomes/${params.genome}/"
params.genome_pattern = "*.{fna,fa}"
params.annotation_pattern = "*.gtf"
params.genomePath = params.genomeDir + params.genome_pattern
params.annotationPath = params.genomeDir + params.annotation_pattern

params.readsDir = "${launchDir}/data/raw_reads/"
params.read_pattern = "*.fastq.gz"
params.reads = params.readsDir + params.read_pattern

params.dev = false
params.number_of_inputs = 2


println """
        L I M S - R N A S E Q   P I P E L I N E
        =======================================
        
        Core Nextflow Options
        launchDir     : ${launchDir}
        workDir       : ${workDir}
        projectDir    : ${projectDir}

        Input Options
        readsDir      : ${params.readsDir}
        genomeDir     : ${params.genomeDir}

        read_pattern  : ${params.read_pattern}
        genome_pattern: ${params.genome_pattern}
        annotation    : ${params.annotationPath}
        """
        .stripIndent()

def get_read_metadata(it) {
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
    .map{get_read_metadata(it)}
    .set{raw_reads}

Channel
    .fromPath(params.genomePath)
    // Takes only the first genome it finds in param.genome
    .first()
    .set{genome}

Channel
    .fromPath(params.annotationPath)
    .first()
    .set{annotation}

process FASTQC {
    tag "${meta.id}"
    // Copies outputs of process to publishdir. No need to use absolute paths
    // in shell/script block
    publishDir "${launchDir}/reports/fastqc", mode: 'copy'

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
    publishDir "${launchDir}/reports/fastqc", mode: 'copy'

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
        path: "${launchDir}/data/trim_reads",
        pattern: '*_trimmed.fq.gz',
        mode: 'copy'
        )

    publishDir(
        path:"${launchDir}/reports/trim_galore",
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

    storeDir "${projectDir}/data/genomes/${params.genome}"

    input:
        path(genome)

    output:
        path('index/*.ht2'), emit: index
        file(genome)

    // TODO: provide core argument dynamically
    shell:
        '''
        mkdir index
        hisat2-build -p 12 !{genome} index/!{params.genome}
        '''
}

process ALIGNMENT {
    tag "${meta.id}"

    publishDir(
        path: "${launchDir}/data/alignments",
        pattern: "*.ba[m,i]",
        mode: 'copy'
    )

    publishDir(
        path: "${launchDir}/reports/hisat2",
        pattern: "*.txt",
        mode: 'copy'
    )

    input:
        tuple val(meta), path(reads)
        path(index)

    output:
        path('*.txt')
        path('*.bam'), emit: alignments
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

process FEATURECOUNTS {
    publishDir "${launchDir}/reports/featureCounts", mode: 'copy'

    input:
        path(alignments)
        path(annotation)

    output:
        path("gene_counts.txt")

    shell:
        '''
        featureCounts \
        --primary \
        -s 2 \
        -T 10 \
        -a !{annotation} \
        -o gene_counts.txt \
        !{alignments}
        '''
}

workflow {
    TRIMMING(raw_reads)
    // // Mixing trimmed and untrimmed reads before fastqc
    all_reads = raw_reads.mix(TRIMMING.out.trim_reads)
    FASTQC(all_reads)
    MULTIQC(FASTQC.out.collect())
    INDEXING(genome)
    ALIGNMENT(TRIMMING.out.trim_reads, INDEXING.out.index)
    FEATURECOUNTS(ALIGNMENT.out.alignments.collect(), annotation)
}
