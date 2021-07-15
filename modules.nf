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
    path fwd
    path rvs

  output:
    path "*_raw.txt" 

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
    tuple val(id),path(unpaired_read)
  output:
    tuple val(id), path("*_U.fastq.gz")

  """
  cat $unpaired_read > ${id}_U.fastq.gz
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
    tuple val(id), path(reads_pair)
    tuple val(id), path(reads_unpair)
    
  output:
    path "*.zip"
    

  """
  ${baseDir}/prog/FastQC/fastqc $reads_pair $reads_unpair --noextract --quiet
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
    path fwd_P
    path rvs_P
    path fwd_U
    path rvs_U

  output:
    path "*cleaned.txt" 
  
  """
  #!/bin/bash
  source ${baseDir}/awk.sh 
  paste <(extract_fastqc -b ${fwd_P}) <(extract_fastqc -o ${rvs_P}) <(extract_fastqc -o ${fwd_U}) <(extract_fastqc -o ${rvs_U}) | body sort -k1 > qc_cleaned.txt
    
  """ 
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
  val mergeUnpair

  output:
  path "*qc_table.txt" 
  script:
  if (mergeUnpair){
  """

  cat ${qc_trim} | sort | grep P.fastq.gz > qc_pair
  awk -F '\\t' 'NR>1{for(i=0;i<2;i++)print}' ${trim_sum} | \
  awk -F '\\t' 'NR==FNR{sur[NR]=\$1;next} {\$2 = \$2 " (" sur[FNR] ")"; print}' OFS='\\t' - qc_pair > qc_pair_nh

  cat ${qc_trim} | sort |grep U.fastq.gz > qc_unpair
  awk -F '\\t' 'NR>1{\$5=\$2+\$3; print}' OFS='\\t' ${trim_sum}| \
  awk -F '\\t' 'NR==FNR{sur[NR]=\$5;dro[NR]=\$4;next} {\$2 = \$2 " (" sur[FNR] ")";print \$0, dro[FNR]}' OFS='\\t' - qc_unpair |sed G > qc_unpair_nh

  cat <(echo "${params.qc_header}") qc_pair_nh > qc_pair.txt
  cat <(echo -e "${params.qc_header}\tDropped") qc_unpair_nh  > qc_unpair.txt
  
  timestamp=\$(date '+%H%M%d%m%y')
  paste ${qc_raw} qc_pair.txt qc_unpair.txt | \
  awk -F '\\t' 'NR>1&&\$NF{\$NF= \$2-\$8-\$14 " (" \$NF ")"}{print}' OFS='\\t' > \${timestamp}_qc_table.txt 
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
      \$NF=\$2+\$6-\$10-\$14-\$18-\$22 " (" (\$2+\$6-\$10-\$14-\$18-\$22)/(\$2+\$6)*100 ")";
      } 
      {print}' OFS='\\t' > qc_table.txt
  
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
  errorStrategy 'ignore'

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

