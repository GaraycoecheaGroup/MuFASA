# StReSSMAn
Short read sequencing Somatic Mutation Analysis and Filtering.
This pipeline will map small sequencing reads to a publicly available reference genome with bwa-mem.
After mapping, small somatic variants are called with Strelka and Mutect2. 
SNV's and INDELs are filtered on quality by several filters. 

# Installation Guide

## Nextflow & Conda
Make sure you have [nextflow](https://www.nextflow.io/docs/latest/index.html) available. 
This pipeline is optimized to be run on a cluster with scheduling software like SLURM available.

The Nextflow pipeline primarily uses conda environments. To install the proper tools and dependencies, use the conda env requirements (.yml) files. All Tools and Packages used in the nextflow workflow are specified in these files. 

## workflow configuration

### nextflow specific
Modify the resources.config configuration file with your run information
SLURM is set as default in the nextflow.config file. If you don't have a SLURM setup available, modify this to the appropriate executor for your setup.

### general configuration / Input
samples.csv, comma seperated file with experiment information (sampleid, fastq.R1, fastq.R2) 
If you have multiple lanes for one sample, add a second row for the same sampleid.
tumor_normal_pairs.txt, comma seperated file with information for the somatic variant calling. On row 1, (tumor/treated) sampleid, tumor, nr
on accompanying second row, same nr: (normal/untreated) sampleid, normal, nr
nr increments per sample pairs (two rows)
filter_parameters.txt. This is the configuration used for [FiNGS](https://pubmed.ncbi.nlm.nih.gov/33602113/) 
This tool is used for filtering SNP's after variant calling

## workflow testing
In the TESTDATA dir is a testset with fastq data to test your installation. 
If your installation is correct, you will be able to run this data set by running:

`./run_nextflow.sh`

### workflow output

|   File                       |   location        |   output                                                      | 
|------------------------------|-------------------|---------------------------------------------------------------|
| *.fastqc.html                | ../FASTQC/        | Quality control checks on raw sequence data                   |
| *_collect_align_metrics.txt  | ../FASTQC/        | Quality metrics for illimina read alignments                  |
| *_collect_wgs_metrics.txt    | ../FASTQC/        | Metrics about coverage and performance of (WGS) mapping       |
| *_samtools_markdup.(bam/bai) | ../MAPPING/       | Read mapping output                                           |
| *.somatic.snvs & indel.vcf   | ../Strelka/       | Strelka somatic SNV's  and InDel variants                     |
| *_Mutect.snvs & indel.vcf    | ../Mutect/        | Mutect2 somatic SNV's and InDel variants                      |
| *_PASS.vcf                   | ../Filtered/      | Strelka and Mutect2 PASS filtered SNV's and InDel variants    |
| *_shared.snvs & indel.vcf    | ../Filtered/      | Intersected results for Mutect/Strelka                        |
| *_pon.snvs & indel.vcf       | ../Filtered/      | Customized filtering for Panel of Normal samples              |
| *.filtered.vcf               | ../snvs_Filtered/ | FiNGS filtered SNV's                                          |
| *.plots.pdf                  | ../snvs_Filtered/ | FiNGS report                                                  |
 


## FASTQ processing
All paired-end reads are assesed on read quality and multiple sequencing runs are merged before trimming 
[cutadapt](https://cutadapt.readthedocs.io/en/stable/) is used to find and remove adapter sequences, primers and poly-A tails. 

## FASTQC 
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Picard CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard) and [Picard CollectAlignmentSummaryMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360040507751-CollectAlignmentSummaryMetrics-Picard) produce reports that can be combined by running:

`multiqc ./FASTQC/`


## MAPPING processing
All reads are mapped against 


## Variant Calling

## Variant Filtering

[^]: Tissue-specific mutagenesis from endogenous guanine damage is suppressed by Polk and DNA repair

