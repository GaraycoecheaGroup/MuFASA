process Strelka {

    label 'strelka'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/strelka'
    publishDir "${params.strelka}"
    //executor 'slurm'

    input:
        tuple val(tumorID), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)

    output:
        val(tumorID) 

    script:
      """
      configureStrelkaSomaticWorkflow.py --normalBam ${normal_bam} --tumorBam ${tumor_bam} --referenceFasta ${params.ref} --outputCallableRegions --runDir ${params.strelka}/${tumorID}_runworkflow/
      python2 ${params.strelka}/${tumorID}_runworkflow/runWorkflow.py -j ${task.cpus} -m local -g 10
      """
}

process handleStrelka {
    label 'mutect_flag'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
    publishDir "${params.strelka}"  
    //executor 'local'

  input:
    val(sampleid)
    
  output:
    tuple val(sampleid), path("*indels.vcf.gz"), path("*snvs.vcf.gz")

  script:
      """
      cp ${params.strelka}/${sampleid}_runworkflow/results/variants/somatic.indels.vcf.gz ${sampleid}.somatic.indels.vcf.gz
      cp ${params.strelka}/${sampleid}_runworkflow/results/variants/somatic.snvs.vcf.gz ${sampleid}.somatic.snvs.vcf.gz
      """
}

process Mutect2 {

  label 'mutect'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/hpc/hub_garayco/software/miniconda3/envs/gatk4'
  //publishDir "${params.mutect}"  
  //executor 'local'

  input:
    //[[a1d2, a1, a1d2_markdups.CRAM, a1_markdups.CRAM, 12]  
    tuple val(tumorid), val(normalid), val(tumor_bam), val(normal_bam), val(chrID)

  output:
    tuple val(tumorid),path("${tumorid}_${chrID}.vcf.gz"), val(chrID)

  script:
      """
      gatk Mutect2 --native-pair-hmm-threads ${task.cpus} -R ${params.ref} --callable-depth ${params.callable_depth} -I ${tumor_bam} -normal ${normalid} -I ${normal_bam} -L chr${chrID} -O ${tumorid}_${chrID}.vcf.gz
      """
}

process Mutect2_flag {

  label 'mutect_flag'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/hpc/hub_garayco/software/miniconda3/envs/gatk4'
  //publishDir "${params.mutect}"  
  //executor 'local'

  input:
    tuple val(tumorid),val(variants),val(chrID)

  output:
    tuple val(tumorid),path("${tumorid}_${chrID}.flagged.vcf.gz"), val(chrID)

  script:
      """
      gatk FilterMutectCalls -R ${params.ref} -V ${variants} -O ${tumorid}_${chrID}.flagged.vcf.gz
      """
}

process Mutect2_concat {
  label 'mutect_flag'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
  publishDir "${params.mutect}"  
  //executor 'local'

  input:
    tuple val(tumorid), path(flagged_files), val(chrIDs)

  output:
    tuple val(tumorid), path("*indels.vcf.gz"), path("*snvs.vcf.gz")

  script:
      """
      bcftools concat -O v -o ${tumorid}_final.vcf ${flagged_files}
      bcftools view -v snps -O z -o ${tumorid}_Mutect.snvs.vcf.gz ${tumorid}_final.vcf
      bcftools view -v indels -O z -o ${tumorid}_Mutect.indels.vcf.gz ${tumorid}_final.vcf
      """
}