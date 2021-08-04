/*
 * Process definitions
 */

/*
 * Step 1a: Read QC before trimming using fastqc
 * Input: Raw reads
 * Output: [samplename].fastqc.zip
 */

process FASTQC_BEFORE_TRIM {
  publishDir "$params.results/fastqc/raw", mode: 'copy'
  tag "$id"

  input:
    tuple val(id), path(reads)

  output:
    path "*.zip"
  
  """
  ${baseDir}/prog/FastQC/fastqc $reads --noextract --quiet
  """
}

/*
 * Step 1b: Extract QC data of each sample 
 * Input: [samplename].fastqc.zip
 * Output: [samplename]_raw.txt
 */

process EXTRACT_FASTQC_BEFORE_TRIM {
  publishDir "$params.results/qc_table", mode: 'copy'
  input: 
    path fastqc_zip
  output:
    path "*_raw.txt" 

  script:
  def fwd = fastqc_zip.findAll { it =~ /.*_1_fastqc.zip/ }.join(' ')
  def rvs = fastqc_zip.findAll { it =~ /.*_2_fastqc.zip/ }.join(' ')
  """
  #!/bin/bash
  source ${baseDir}/awk.sh 
  paste <(extract_fastqc -b ${fwd}) <(extract_fastqc -o ${rvs}) | body sort -k1 > qc_raw.txt
  """
}

/*
 * Step 1c: Merge all sample's QC data into one table
 * Input: [samplename]_raw.txt
 * Output: qc_table_raw.txt
 */

process CREATE_QCTABLE_BEFORE_TRIM {
  publishDir "$params.results/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "qc_raw.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "${params.qc_header}") qc_table > qc_raw.txt
  """


}
/*
 * Step 1d: Summarise fastqc results using multiqc
 * Input: [samplename].fastqc.zip
 * Output: "raw_data/ , raw.html"
 */
process MULTIQC_FASTQC_BEFORE_TRIM {
  publishDir "$params.results/multiqc", mode: 'copy'
  errorStrategy 'ignore'

  input:
    path fastqc_before_all

  output:
    path "raw.html"
    path "raw_data"

  """
  multiqc ${fastqc_before_all} --interactive --filename raw
  """
}


/*
 * Step 2a. Trim adapters and low quality reads
 * Input: Raw reads
 * Output: "[samplename]_1P.fq.gz, [samplename]_1U.fq.gz, [samplename]_2P.fq.gz, [samplename]_2U.fq.gz, [samplename]_summary.txt" 
 */

process TRIM {
  publishDir "$params.results/trim_results", mode: 'copy'
  tag "$id"
 // errorStrategy 'ignore'

  input:
    tuple val(id), path(reads)
    val qc_option

  output:
    tuple val(id), path("*P.fastq.gz")
    tuple val(id), path("*U.fastq.gz")
    path "*_summary.txt"
  shell:
  '''
  echo !{reads[0]} !{reads[1]}
  java -jar !{baseDir}/prog/trimmomatic/trimmomatic-0.39.jar PE -phred33 -trimlog !{id}_log.txt -summary !{id}_summary.txt \
  -basein !{reads[0]} \
  -baseout !{id}.fastq.gz \
  ILLUMINACLIP:!{baseDir}/prog/trimmomatic/adapters/adapter.fa:2:30:10 !{qc_option} 
  '''
}

/*
 * Step 2b: Extract cleaning result from each sample
 * Input: [samplename]_summary.txt
 * Output: [samplename]_ext.txt
 */
process EXTRACT_TRIM_LOG {
  input:
  path summarylog

  output:
  path "*_ext.txt"

  """
  cat ${summarylog} | awk -F ':' '{print \$2}' | paste -s | 
  awk  '{print \$3,\$5,\$7,\$9}' OFS='\\t' > ${summarylog.baseName}_ext.txt
  """
}

/*
 * Step 2c: Merge cleaning result into one table
 * Input: [samplename]_ext.txt
 * Output: trim_log.txt
 */

process CREATE_TRIM_SUMMARY_TABLE {
  publishDir "$params.results/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "trim_summary.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "${params.trim_header}") qc_table > trim_summary.txt

  """
}

process MERGE_UNPAIRED_READ {
 
  input: 
    tuple val(id),path (unpaired_read)
  output:
    tuple val(id),path("*_U.fastq.gz")

  script: 
    def unpair = unpaired_read.findAll{ it =~ /.*U.fastq.gz/}.join(' ')
      """
      cat $unpair > ${id}_U.fastq.gz
      """
}

/*
 * Step 3a: Read QC after trimming using fastqc
 * Input: [samplename]_{1,2}P.fq.gz, [samplename]_{1,2}U.fq.gz
 * Output: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 */
process FASTQC_AFTER_TRIM {
  publishDir "$params.results/fastqc/trimmed", mode: 'copy'
  tag "$id"

  input:
    tuple val(id), path(pair)
    tuple val(id_u), path(unpair)
  output:
    path "*.zip"
    

  """
  ${baseDir}/prog/FastQC/fastqc $pair $unpair  --noextract --quiet
  """
}


/*
 * Step 3b: Extract QC data from cleaned samples 
 * Input: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 * Output: [samplename]_{1,2}P_Trimmed.txt , [samplename]_{1,2}U_Trimmed.txt
 */

process EXTRACT_FASTQC_AFTER_TRIM {
  publishDir "$params.results/qc_table", mode: 'copy'
  input: 
    path fastqc_zip
  output:
    path "*cleaned.txt" 
  
  script:
    def fwd_pair = fastqc_zip.findAll { it =~ /.*1P_fastqc.zip$/ }.join(' ')
    def rvs_pair = fastqc_zip.findAll { it =~ /.*2P_fastqc.zip$/ }.join(' ')
    if (params.mergeUnpair){
    def unpair = fastqc_zip.findAll {it =~ /.*U_fastqc.zip$/ }.join(' ')
      """
      #!/bin/bash
      source ${baseDir}/awk.sh 
      paste <(extract_fastqc -b ${fwd_pair}) <(extract_fastqc -o ${rvs_pair}) <(extract_fastqc -o ${unpair}) | body sort -k1 > qc_cleaned.txt
        
      """
    } else {
    def fwd_unpair = fastqc_zip.findAll { it =~ /.*1U_fastqc.zip$/ }.join(' ')
    def rvs_unpair = fastqc_zip.findAll { it =~ /.*2U_fastqc.zip$/ }.join(' ')
      """
      #!/bin/bash
      source ${baseDir}/awk.sh 
      paste <(extract_fastqc -b ${fwd_pair}) <(extract_fastqc -o ${rvs_pair}) <(extract_fastqc -o ${fwd_unpair}) <(extract_fastqc -o ${rvs_unpair}) | body sort -k1 > qc_cleaned.txt
        
      """ 
  
    }
}
/*
 * Step 3c: Merge QC data of paired and unpaired samples into two different table
 * Input: [samplename]_{1,2}P_Trimmed.txt , [samplename]_{1,2}U_Trimmed.txt
 * Output: "qc_pair.txt, qc_unpair.txt"
 */

process FINALISE_QCTABLE {
  publishDir "$params.results/qc_table", mode: 'copy'
  input:
  path qc_raw
  path qc_trim
  path trim_sum

  output:
  path "*qc_table.txt" 
  script:
  if (params.mergeUnpair){
  """
  #!/bin/bash 
  join --header -t \$'\\t' qc_raw.txt qc_cleaned.txt | \
  awk -F '\t' 'NR==1{
      \$(NF+1)="Drop";
      } 
  NR>1{
      \$10= \$10 " (" \$10/\$2*100 ")";
      \$14= \$14 " (" \$14/\$2*100 ")"; 
      \$18= \$18 " (" \$18/\$2*100 ")"; 
      \$(NF+1)=\$2+\$6-\$10-\$14-\$18-\$22 " (" (\$2+\$6-\$10-\$14-\$18)/(\$2+\$6)*100 ")";
      } 
      {print}' OFS='\\t' > qc_table.txt
  
  """
 
  } else {
  
  """
  #!/bin/bash 
  join --header -t \$'\\t' qc_raw.txt qc_cleaned.txt | \
  awk -F '\t' 'NR==1{
      \$(NF+1)="Drop";
      } 
  NR>1{
      \$10= \$10 " (" \$10/\$2*100 ")";
      \$14= \$14 " (" \$14/\$2*100 ")"; 
      \$18= \$18 " (" \$18/\$2*100 ")"; 
      \$22= \$22 " (" \$22/\$2*100 ")"; 
      \$(NF+1)=\$2+\$6-\$10-\$14-\$18-\$22 " (" (\$2+\$6-\$10-\$14-\$18-\$22)/(\$2+\$6)*100 ")";
      } 
      {print}' OFS='\\t' > qc_table_tmp.txt
  echo -e "RAW\t\t\t\t\t\t\t\t\tTRIMMED" > header
  echo -e "FWD\t\t\t\t\tRVS\t\t\t\tPAIR_FWD\t\t\t\tPAIR_RVS\t\t\t\tUNPAIR_FWD\t\t\t\tUNPAIR_RVS" >> header
  cat header qc_table_tmp.txt > qc_table.txt
  """
  }
}

/*
 * Step 3d: Summarise fastqc results using multiqc
 * Input: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 * Output: "cleaned_data/ , cleaned.html"
 */

process MULTIQC_FASTQC_AFTER_TRIM {
  publishDir "$params.results/multiqc", mode: 'copy'
  // errorStrategy 'ignore'

  input:
    path fastqc_after_all

  output:
    path "trimmed.html"
    path "trimmed_data"

  """
  multiqc ${fastqc_after_all} --interactive --filename trimmed 
  """
}

/*
 * Step 4: Merge extracted QC data from both before and after cleaning into one final table
 * Input: qc_table_raw.txt, qc_pair.txt, qc_unpair.txt
 * Output: qc_table.txt
 */

process MERGE_QCTABLE {
  publishDir "$params.results/qc_table", mode: 'copy'

  input:
    path raw_qc_table
    path trim_qc_table_pair
    path trim_qc_table_unpair

  output:
    path "*qc_table.txt"
  
  """ 
  timestamp=\$(date '+%d%m%y_%H%M%S')
  paste ${raw_qc_table} ${trim_qc_table_pair} ${trim_qc_table_unpair} > \${timestamp}_qc_table.txt
  """
}

