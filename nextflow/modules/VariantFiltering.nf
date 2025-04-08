process PASSfilter {
  label 'merge_label'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
  publishDir "${params.filtered_dir}"  

  input:
    tuple(val(sample_id), path(raw_vcf), val(tool), val(type)) //-> change to path..

  output:
    tuple(val(sample_id),path("*_PASS.vcf.g*"),val(tool), val(type))

  script:
      """
      bcftools view --threads ${task.cpus} -O z -o ${sample_id}.${tool}.${type}_PASS.vcf.gz -f PASS ${raw_vcf}
      bcftools index ${sample_id}.${tool}.${type}_PASS.vcf.gz 
      """

}

process ISEC_overlap {
    label 'isec_overlap'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
    publishDir "${params.filtered_dir}"  
    //executor 'local'

    input:
      tuple(val(sample_id),path(strelka),path(strelka_index),path(mutect), path(mutect_index), val(type))

    output:
      tuple(val(sample_id),path("${sample_id}_shared_${type}.vcf"), val(type))

    script:
      """
      bcftools isec --threads ${task.cpus} -p isec_${sample_id} ${mutect} ${strelka}
      cp isec_${sample_id}/0002.vcf ${sample_id}_shared_${type}.vcf
      """
}

process PoN_filter {
    label 'isec_overlap'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
    publishDir "${params.filtered_dir}" 

    input:
      tuple(val(sample_id),path(shared),val(type))

    output:
      tuple(val(sample_id),path("${sample_id}_pon_${type}.vcf"), val(type))

    script:
      """
      if [ ! -f ${shared}.gz.csi ]; then
        bgzip ${shared}
        bcftools index ${shared}.gz
      fi

      bcftools isec --threads ${task.cpus} -O v -p isec_${sample_id} -C ${shared}.gz ${params.pon}
      cp isec_${sample_id}/0000.vcf ${sample_id}_pon_${type}.vcf

      """
}

process FiNGS {
    label 'fings'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/fings'
    publishDir "${params.snvs_filtered}" 

    input:
         tuple(val(sample_id),path(pon),path(tumorbam),path(tumorbai),path(normalbam), path(normalbai))

    output:
        val(sample_id)

    script:
      """
        fings -j ${task.cpus} --overwrite -r ${params.ref} -p ${params.fings} -n ${normalbam} -t ${tumorbam} -v ${pon} -d ${params.snvs_filtered}/${sample_id}_fings/ --PASSonlyin --PASSonlyout

      """
}

process handleFings {
    label 'plot_vaf'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/envs/pipeline'
    publishDir "${params.snvs_filtered}"  
    //executor 'local'

    input:
        val(sample_id)
      
    output:
        val(sample_id)

    script:
      """ 
        cp ${params.snvs_filtered}/${sample_id}_fings/*.filtered.vcf ${params.snvs_filtered}/${sample_id}_fings.vcf
        cp ${params.snvs_filtered}/${sample_id}_fings/plots.pdf ${params.pdf_path}/${sample_id}_fings.pdf
      """
}

process snvsVAF {
    label 'vaf'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/pipeline'
    publishDir "${params.snvs_filtered}"  
    //executor 'local'

    input:
         val(sample_id)

    output:
        tuple(val(sample_id),path("${sample_id}.vaf.vcf")) 

    script:
      """
      python ${params.script_dir}/filter_mutect2_allelic_frequency.py ${params.snvs_filtered}/${sample_id}_fings.vcf ${sample_id}.ALLvaf.tab ${sample_id}.vaf.vcf ${sample_id}.PASSvaf.tab ${params.minvaf} ${params.maxvaf}

      """
}

process indelFiltering {
    label 'indel'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/hpc/hub_garayco/software/miniconda3/pipeline'
    publishDir "${params.indel_filtered}" 

    input:
        tuple val(sample_id),path(pon),val(type) 

    output:
        tuple(val(sample_id),path("${sample_id}.NoAlt.vcf"),path("${sample_id}.SBtolerant.vcf"), path("${sample_id}.SBzero.vcf"), path("${sample_id}.RF.vcf"),path("${sample_id}.AF.vcf"),path("${sample_id}.AF.csv"))

    script:
      """
      python ${params.script_dir}/filter_shared_indels.py ${pon} ${sample_id} ${params.minvaf} ${params.maxvaf}

      """

}






