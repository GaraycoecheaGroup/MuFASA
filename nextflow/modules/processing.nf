
def defineOrderISEC(grouped){
    reordered = grouped
        .map{ it ->
            if (it[2][0] == "mutect"){
                //println "mutect, strelka"
                [it[0],it[1][1][0],it[1][1][1],it[1][0][0],it[1][0][1],it[3]] 
            }else{
               //println "strelka, mutect"
                [it[0],it[1][0][0],it[1][0][1],it[1][1][0],it[1][1][1],it[3]]
            }//.view()
        }
    reordered
}



def preVariantCalling(paired_channel) {
    tumor_normals = paired_channel
        .map{ item -> 
                if (item[3][0] == "normal"){
                    //println "wrong order"
                    [[item[0][1],item[0][0]],[item[1][1],item[1][0]],[item[2][1],item[2][0]],[item[3][1],item[3][0]],item[4]]
                }else{
                    [item[0],item[1],item[2],item[3],item[4]]
                } 
            }//.view()
    tumor_normals
}

def add_indel_snvs_tags(transposed){
    added_tag = transposed
        .map { item ->
                if (item[1] =~ "somatic.indels.vcf.gz"){
                    [item[0],item[1],"strelka","indel"]
                }else if (item[1] =~ "somatic.snvs.vcf.gz"){
                    [item[0],item[1],"strelka", "snvs"]
                }else if (item[1] =~ "Mutect.indels.vcf.gz"){
                    [item[0],item[1],"mutect", "indel"]
                }else if (item[1] =~ "Mutect.snvs.vcf.gz"){
                    [item[0],item[1],"mutect", "snvs"]
                }
            }
    added_tag
}



def getSampleIds(dir) {
    dir = dir.tokenize().collect{"$it/*_1.fq.gz"}   
    Channel
    .fromPath(dir, type:'file')
    .ifEmpty { error "No R1 fastq.gz files found in ${dir}." }
    .map { r1_path ->
        sample_id = r1_path.getSimpleName().split('_')[0]
    }
}



def flowcellLaneFromFastq(path) {                                                                                                                                      
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab                                                           
                                                                                                                                                                       
    // parse first line of a FASTQ file (optionally gzip-compressed)                                                                                                   
    // and return the flowcell id and lane number.                                                                                                                     
    // expected format:                                                                                                                                                
    // xx:yy:FLOWCELLID:LANE:... (seven fields)                                                                                                                        
    // or                                                                                                                                                              
    // FLOWCELLID:LANE:xx:... (five fields)                                                                                                                            
    InputStream fileStream = new FileInputStream(path.toFile())                                                                                                        
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)                                                                                             
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')                                                                                                        
    BufferedReader buffered = new BufferedReader(decoder)                                                                                                              
    def line = buffered.readLine()                                                                                                                                     
    assert line.startsWith('@')                                                                                                                                        
    line = line.substring(1)                                                                                                                                           
    def fields = line.split(' ')[0].split(':')                                                                                                                                                                                                                                                                        
    String fcid                                                                                                                                                        
    int lane                                                                                                                                                           
                                                                                                                                                                       
    if (fields.size() == 7 || fields.size() == 8) {                                                                                                                    
        // CASAVA 1.8+ format                                                                                                                                          
        fcid = fields[2]                                                                                                                                               
        lane = fields[3].toInteger()                                                                                                                                   
    }                                                                                                                                                                  
    else if (fields.size() == 5) {                                                                                                                                     
        fcid = fields[0]                                                                                                                                               
        lane = fields[1].toInteger()                                                                                                                                   
    }                                                                                                                                                                 
    [fcid, lane]                                                                                                                                      
}                      

def extractAllFastqFromDir(dir) {
    // Original code from: https://github.com/SciLifeLab/Sarek - MIT License - Copyright (c) 2016 SciLifeLab
    dir = dir.tokenize().collect{"$it/**_R1.fastq.gz"}
    Channel
    .fromPath(dir, type:'file')
    .ifEmpty { error "No R1 fastq.gz files found in ${dir}." }
    .map { r1_path ->
        fastq_files = [r1_path]
        sample_id = r1_path.getSimpleName().split('_')[0]
        r2_path = file(r1_path.toString().replace('_R1', '_R2'))
        if (r2_path.exists()) fastq_files.add(r2_path)
        (flowcell, lane) = flowcellLaneFromFastq(r1_path)
        rg_id = "${sample_id}_${flowcell}_${lane}"
        [sample_id, rg_id,fastq_files]
    }
}
