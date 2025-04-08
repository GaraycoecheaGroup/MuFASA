process WgsMetrics {
    label 'wgs_metrics'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/gatk4'
    publishDir params.fastqc_path, mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple(val(sample_id), path(input))

    output:
        tuple val(sample_id), path("${sample_id}_collect_wgs_metrics.txt")


    script:
        """
        picard CollectWgsMetrics R=${params.ref} I=${input} O=${sample_id}_collect_wgs_metrics.txt VALIDATION_STRINGENCY=SILENT
    
        """
}

process AlignmentSummary {
    label 'Align_sum'
    shell = ['/bin/bash', '-euo', 'pipefail']    
    conda '/hpc/hub_garayco/software/miniconda3/envs/gatk4'
    publishDir params.fastqc_path, mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple(val(sample_id), path(input))

    output:
         tuple(val(sample_id), path("${sample_id}_collect_align_metrics.txt"))

    script:
        """
        picard CollectAlignmentSummaryMetrics R=${params.ref} I=${input} O=${sample_id}_collect_align_metrics.txt VALIDATION_STRINGENCY=SILENT
        
        """

}

process PlotVAF {
    label 'plot_vaf'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
    errorStrategy 'ignore'
    publishDir "${params.pdf_path}"

    input:
        tuple(val(sample_id), path(input), val(type), val(color))

    output:
        tuple(val(sample_id), path("${sample_id}_${type}_VAFplot.pdf"))

    script:
        """
        python ${params.script_dir}/plot_VAF_mutect.vcf.py ${input} ${color} ${type} ${sample_id}
        """
}