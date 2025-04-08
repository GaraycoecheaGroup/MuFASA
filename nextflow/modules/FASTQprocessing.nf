process MergeFastq {

    label 'merge_label'
    conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
    publishDir params.merged_dir

    input:
        tuple val(sample_id), path(fastq_files), val(pair)
    
    output:
        tuple(val(sample_id), path("*fq.gz")) 

    script:
    """
    cat ${fastq_files} > ${sample_id}_${pair}.fq.gz 

    """
}


process FastQC {
    label 'FastQC_label'
	conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
    errorStrategy 'ignore'

    input:
        tuple val(sample_id), path(fastq)

    script:
    """
    fastqc -t ${task.cpus} -o ${params.fastqc_path} ${fastq}

    """
}

process CUTadapt {
    label 'cutadapt'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'

    input:
        tuple(val(sample_id), path(fastq))

    output:
        tuple(val(sample_id), path("trimmed_*fq.gz"))

    script:
        """
        cutadapt -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" -o trimmed_${sample_id}_1.fq.gz -p trimmed_${sample_id}_2.fq.gz --max-n 15 --max-ee 24 --cores=${task.cpus} ${fastq}
        """
}

