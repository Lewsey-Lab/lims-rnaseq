# lims-rnaseq

## Introduction
lims-rnaseq is a simple rnaseq pipeline adapted to the older linux kernel on
the La Trobe Institute for Molecular Science (LIMS) high performance computing
cluster (HPCC). The older linux kernel is unable to run the mainstream
[nf-core/rnaseq](https://nf-co.re/rnaseq) due to its usage of more recent linux
features/programs.

As such, this pipeline was written to be compatible with what is currently
available on the LIMS HPCC.

## Pipeline Summary
1. Adapter and quality trimming (`trimgalore/0.6.3`)
2. Read QC (`fastqc/0.11.9`, `multiqc/1.9`)
3. Creation of genome index (`hisat-gcc/2.1.0`)
4. Read alignment and indexing (`hisat-gcc/2.1.0`, `samtools-gcc7/1.9`)
5. Gene-level quantification (`featureCounts 1.6.0`, `subread-gcc/1.6.0`)

## Installation
1. Install [java 8 or later (up to
   15)](https://www.oracle.com/java/technologies/javase-downloads.html) if
   installing nextflow on your home PC. If working on the HPCC, load the java
   1.8 module:
   ```
   module load java/1.8.0_66
   ```
2. Install
   [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)
   (`>=21.04.0`) using the following commands:
   ```
   wget -qO- https://get.nextflow.io | bash
   ```
3. Move the newly created executable file `nextflow` to a folder accessible by
   `PATH`.
   
4. Pull the `lims-rnaseq` pipeline with the following command: 
   ```
   nextflow pull https://github.com/SpikyClip/lims-rnaseq
   ```
   This stores the pipeline in `.nextflow/assets/SpikyClip/lims-rnaseq`
   allowing it to be used in any directory.

## Quick Start

1. Create a new folder for your project and `cd` into it:
   ```
   mkdir <my_project>
   cd <my_project>
   ```
2. Create the following folder structure for inputs:
   ```
   data
   ├── genomes
   │   └── <genome_name>
   │       ├── <genome_fasta_or_fna>.fna
   │       └── <annotation_file>.gtf
   └── raw_reads
       ├── <accession_id>_<this_is_a_SE_read>.fastq.gz
       ├── <accession_id>_<this_is_a_PE_read>_1.fastq.gz
       └── <accession_id>_<this_is_a_PE_read>_2.fastq.gz
   ```
   Take note of the name of the genome folder, `<genome_name>`. This will be
   later used as an argument to direct nextflow to the `.fna/.fa` and `.gtf`
   files (there should only be one `.fna/.fa file` and one `.gtf` file per
   genome folder as this pipeline sticks to one genome per run). Also take note
   of the format of the read names. PE reads are automatically distinguished by
   an identical name apart from the trailing `_1` and `_2`. If the body name is
   not identical, it would be treated as a SE read. (**NOTE:** I have yet to
   write a test to make sure inputs are correct, but will do so in the future.
   **All reads are also treated as reverse reads**, will probably add a way to
   detect it from the file name.)

3. Load the java module (if not already loaded) and run the pipeline in the
   `<my_project>` directory:
   ```
   module load java/1.8.0_66

   nextflow \
   run lims-rnaseq \
   -profile cluster \
   --genome <genome_name> \
   ```
   Omit the `-profile` argument if not running on the LIMS-HPCC.
4. Stop the pipeline at any time with ctrl + c without losing progress. To
   resume the last run, add the optional argument `-resume`. Results get output
   into the `data` folder and diagnostics into the `reports` folder

## How does it work?

### Results
Results for the user are published to the project folder under newly created
subfolders in `data` and `reports`.

### `work` folder
Nextflow creates a `work` folder containing checksummed folders where the
actually process files are produced. These folders are not meant to be touched
manually (but can be useful for debugging). These folders allow the pipeline to
recall its process, allowing it to stop and resume without losing information.
It also ensures that there will not be any cross process collisions as each
folder is unique. It relies on checksums from inputs and outputs to detect if
anything has changed. To learn more, check out this
[post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

### Genome index files

Genome indexes have been coded to be stored in a more permanent place in the
original project folder in
`.nextflow/assets/SpikyClip/data/genomes/<genome_name>` in the user directory.
This allows you to use the same stored indexes without re-running the indexes,
even in another project directory. **NOTE: You still need the same
`<genome_fasta_or_fna>.fna` in your project `data/genomes/<genome_name>` folder
as this is how nextflow links the input to the correct indexes** (I'm working
on a solution that only relies on a `--genome <genome_name>` argument)

### Processes

Processes are done in parallel whenever possible (e.g. fastqc and trimming on
multiple files). The `-profile cluster` argument tells nextflow to use the
module loading system on the HPCC to load the necessary modules for each
process, and also to generate a slurm job for each process. If you plan on
running a long nextflow job on the HPCC, wrap the `nextflow run` command in a
shell script and submit it like any other slurm job.

### Debugging

If you run into errors, checkout the `.nextflow.log` file created in the
project folder for debugging. Once the errors are fix, run the pipeline with
the `-resume` to continue work. As the pipeline is a work in progress, please
report any issues or recommendations
[here](https://github.com/SpikyClip/lims-rnaseq/issues).

### Pipeline updates

To update the pipeline run:
```
nextflow pull lims-rnaseq
```

## TODO
1. Allow specification of stored genome index without requiring .fna/.fa file.
2. Optimise CPU/Memory usage for each process.
3. Get email notifications working.
4. Automatic removal of empty gene_id in .gtf before featureCounts?
5. Split up processes into sub-workflows

## Project Organization

Here is a diagram of the pipeline package structure in
`.nextflow/assets/SpikyClip/lims-rnaseq`

    lims-rnaseq
    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   └── genomes        <- Folder containing stored genome index folders
    ├── main.nf            <- Main script executed by Nextflow
    ├── nextflow.config    <- Nextflow configuration file
    ├── pipeline_info      <- Folder containing pipeline info and diagrams
    |
    └── src                <- Source code for use in this project
        ├── data           <- TODO: Scripts to download or generate data
        │
        ├── features       <- TODO: Scripts to turn raw data into features for modeling
        |
        ├── models         <- TODO: Scripts to train models to make predictions
        |
        └── visualization  <- TODO: Scripts to create exploratory and results oriented visualizations
