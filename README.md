# Cut&Run Nextflow Pipeline
Cut&Run;Run alignment and bed fragment bed generation pipeline
Cut&Run Nextflow Pipeline

This repository hosts a Nextflow DSL2 pipeline for processing Cut&Run sequencing data from raw FASTQ reads through trimmed, aligned, and fragment-length analysis stages. The workflow is designed for HPC or cloud environments and leverages your cluster’s module system.

Table of Contents

Features

Requirements

Installation

Configure

Usage

Input & Output

Pipeline Outline

Process Descriptions

Profiles

Example

Troubleshooting

License

Features

Adapter trimming via Cutadapt

Alignment to reference genome with Bowtie2

Post-alignment processing with Samtools & Bedtools

Fragment conversion (BEDPE to BED)

Fragment-length distribution computation

Compression & de-duplication steps

Module-based for HPC environments (SLURM, LSF)

Requirements

Nextflow (>=24.x)

Unix-like environment

HPC modules: cutadapt, bowtie2, samtools, bedtools, perl

Slurm or LSF scheduler (profiles provided)

Installation

Clone this repo:

git clone https://github.com/yourorg/cutr_pipeline.git
cd cutr_pipeline

Ensure Nextflow is installed and in your $PATH:

curl -s https://get.nextflow.io | bash
mv nextflow ~/bin/

Configure

Edit nextflow.config to set:

params.DATA: path to input FASTQ directory

params.RESULTS: output directory

params.SAMPLES: list of sample IDs

params.GENOME: Bowtie2 index prefix

params.SCRIPTS: path to helper scripts

Resource directives (memory, time, queue) under your chosen profile

Usage

Run the pipeline with your chosen scheduler profile. For example, on a Slurm cluster:

nextflow run main.nf -c nextflow.config -profile slurm

Input & Output

Input: Paired FASTQ files named as ${sample}_R1_*.fastq.gz and ${sample}_R2_*.fastq.gz in params.DATA.

Output: For each sample, a subdirectory under params.RESULTS containing:

Trimmed FASTQ files

Alignment SAM/BAM/BED files

Fragment-length distributions (*_braw.len, *_bnorm.len)

Compressed BED files and deduplicated outputs

Pipeline Outline

samples_ch -> fastqs_ch
fastqs_ch → cutadapt → trimmed_ch
trimmed_ch → bowtie2 → aligned_ch
aligned_ch → postAlign → pebeds_ch
pebeds_ch → lengthFiles → lengths_ch
pebeds_ch → zipBeds → gzbed_ch
gzbed_ch → uniqBeds → unique_ch

Process Descriptions

Process

Description

cutadapt

Trims adapters and filters reads by length

bowtie2

Aligns trimmed reads to the genome, outputs SAM

postAlign

Sorts SAM to BAM, filters properly-paired reads, converts to BEDPE and then BED

pebeds

Converts BEDPE to single-end BED fragments via Perl helper (bedpe2bed.pl)

lengthFiles

Computes raw counts and normalized fragment-length distributions

zipBeds

Compresses fragment BED files with gzip

uniqBeds

Removes duplicate fragments via uniq_STDOUT.pl and outputs gzipped BED

Profiles

local: run everything locally

slurm: submit via Slurm (--executor slurm, set process.queue, process.time, etc.)

lsf: submit via LSF (--executor lsf)

Example

Prepare your config:

params {
  DATA    = 'data/'          
  RESULTS = 'results/'       
  SAMPLES = ['sample1','sample2']
  GENOME  = 'hg38/index'     
  SCRIPTS = 'code/'          
}
profiles {
  slurm { process.queue='compute'; process.time='48h' }
}

Launch:

nextflow run main.nf -c nextflow.config -profile slurm

Troubleshooting

invalid partition: verify process.queue matches your cluster’s partitions (see sinfo).

QOS errors: set process.memory and clusterOptions (--qos) to satisfy your policies.

Syntax errors: ensure DSL2 is enabled (nextflow.enable.dsl = 2) and braces match.

License

This pipeline is released under the MIT License.

Maintained by kiranvpaul
