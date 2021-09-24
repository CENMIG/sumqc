/*
 * Process definitions
 */

/*
 * Read QC before cleaning using fastqc
 * Input: Raw reads
 * Output: [samplename].fastqc.zip
 */

process FASTQC_BEFORE_CLEANING_SE {
  publishDir "$params.output/fastqc/raw", mode: 'copy'
  tag "$reads"

  input:
    path reads

  output:
    path "*.zip"
  
  """
  ${baseDir}/prog/FastQC/fastqc $reads --noextract --quiet
  """
}

/*
 * Read QC before cleaning using fastqc
 * Input: Raw reads
 * Output: [samplename].fastqc.zip
 */

process FASTQC_BEFORE_CLEANING_PE {
  publishDir "$params.output/fastqc/raw", mode: 'copy'
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
 * Extract FastQC results
 * Input: [samplename]_fastqc.zip
 * Output: [samplename]_raw.txt
 */

process EXTRACT_FASTQC_BEFORE_CLEANING {
  publishDir "$params.output/qc_table", mode: 'copy'
  input: 
    path fastqc_zip
  output:
    path "*_raw.txt" 

  script:
  def fwd = fastqc_zip.findAll { it =~ /.*_1_fastqc.zip/ }.join(' ')
  def rvs = fastqc_zip.findAll { it =~ /.*_2_fastqc.zip/ }.join(' ')
  if (params.SE)
  """
  #!/bin/bash
  source ${baseDir}/awk.sh
  extract_fastqc -b -h ${fastqc_zip} | body sort -k1 > qc_raw.txt
  """
  else if (!params.SE) 
  """
  #!/bin/bash
  source ${baseDir}/awk.sh 
  paste <(extract_fastqc -b ${fwd}) <(extract_fastqc -o ${rvs}) | body sort -k1 > qc_raw.txt
  """
 
}



/*
 * Summarise FastQC results before cleaning using MultiQC
 * Input: [samplename].fastqc.zip
 * Output: "raw_data/ , raw.html"
 */

process MULTIQC_BEFORE_CLEANING {
  publishDir "$params.output/multiqc", mode: 'copy'
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
 * Trim adapters and drop low quality reads
 * Input: Raw reads
 * Output: "[samplename].fq.gz
 */

process CLEAN_SE {
  publishDir "$params.output/cleaned_results", mode: 'copy'
  tag "$reads"
 // errorStrategy 'ignore'

  input:
    path reads
    val qc_options

  output:
    path "*_cleaned.fastq.gz"
    path "*_summary.txt"
    
  when: 
    params.SE

  shell:
  '''
  echo !{reads}
  java -jar !{baseDir}/prog/trimmomatic/trimmomatic-0.39.jar SE -phred33 -trimlog !{reads}_log.txt -summary !{reads}_summary.txt \
  !{reads} \
  !{reads}_cleaned.fastq.gz \
  ILLUMINACLIP:!{baseDir}/prog/trimmomatic/adapters/adapter.fa:2:30:10 !{qc_options} 
  '''
}

/*
 * Trim adapters and drop low quality reads
 * Input: Raw reads
 * Output: "[samplename]_1P.fq.gz, [samplename]_1U.fq.gz, [samplename]_2P.fq.gz, [samplename]_2U.fq.gz, [samplename]_summary.txt" 
 */

process CLEAN_PE {
  publishDir "$params.output/cleaned_results", mode: 'copy'
  tag "$id"
 // errorStrategy 'ignore'

  input:
    tuple val(id), path(reads)
    val qc_options

  output:
    tuple val(id), path("*P.fastq.gz")
    tuple val(id), path("*U.fastq.gz")
    path "*_summary.txt"

  when:
    !params.SE

  shell:
  '''
  echo !{reads[0]} !{reads[1]}
  java -jar !{baseDir}/prog/trimmomatic/trimmomatic-0.39.jar PE -phred33 -trimlog !{id}_log.txt -summary !{id}_summary.txt \
  !{reads[0]} !{reads[1]} \
  -baseout !{id}.fastq.gz \
  ILLUMINACLIP:!{baseDir}/prog/trimmomatic/adapters/adapter.fa:2:30:10 !{qc_options} 
  '''
}

/*
 * Extract read statistics from Trimmomatic log files
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

/* Create read summary statistics table
 * Input: trim_summary.txt
 * Output: trim_summary.txt
 */


process CREATE_TRIM_SUMMARY_TABLE {
  publishDir "$params.output/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "trim_summary.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "BothSurvied\tForwardOnlySurvived\tReverseOnlySurvied\tDroppedRead") qc_table > trim_summary.txt

  """
}

/* Merge unpaired F and R reads' stats 
 * Input: [samplename]_ext.txt
 * Output: trim_log.txt
 */

process MERGE_UNPAIRED_READS {
 
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
 * SE read QC after cleaning using FastQC
 * Input: [samplename]_{1,2}P.fq.gz, [samplename]_{1,2}U.fq.gz
 * Output: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 */
process FASTQC_AFTER_CLEANING_SE {
  publishDir "$params.output/fastqc/cleaned", mode: 'copy'
  tag "$id"

  input:
    path cleaned_reads
  output:
    path "*.zip"
    

  """
  ${baseDir}/prog/FastQC/fastqc $cleaned_reads   --noextract --quiet
  """
}

/*
 * PE read QC after cleaning using fastqc
 * Input: [samplename]_{1,2}P.fq.gz, [samplename]_{1,2}U.fq.gz
 * Output: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 */
process FASTQC_AFTER_CLEANING_PE {
  publishDir "$params.output/fastqc/cleaned", mode: 'copy'
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
 * Extract FastQC results after cleaning
 * Input: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 * Output: [samplename]_{1,2}P_Trimmed.txt , [samplename]_{1,2}U_Trimmed.txt
 */

process EXTRACT_FASTQC_AFTER_CLEANING {
  publishDir "$params.output/qc_table", mode: 'copy'
  input: 
    path fastqc_zip
  output:
    path "*cleaned.txt" 
  
  script:
    def fwd_pair = fastqc_zip.findAll { it =~ /.*1P_fastqc.zip$/ }.join(' ')
    def rvs_pair = fastqc_zip.findAll { it =~ /.*2P_fastqc.zip$/ }.join(' ')
    def fwd_unpair = fastqc_zip.findAll { it =~ /.*1U_fastqc.zip$/ }.join(' ')
    def rvs_unpair = fastqc_zip.findAll { it =~ /.*2U_fastqc.zip$/ }.join(' ')
  if (params.SE){
      """
      #!/bin/bash
      source ${baseDir}/awk.sh
      extract_fastqc -oh ${fastqc_zip} | body sort -k1 > qc_cleaned.txt
      """
  } else if (!params.SE) { 
    if (params.mergeUnpaired){
    def unpair = fastqc_zip.findAll {it =~ /.*U_fastqc.zip$/ }.join(' ')
      """
      #!/bin/bash
      source ${baseDir}/awk.sh 
      paste <(extract_fastqc -b ${fwd_pair}) <(extract_fastqc -o ${rvs_pair}) <(extract_fastqc -o ${unpair}) | body sort -k1 > qc_cleaned.txt
        
      """
    } else {
      """
      #!/bin/bash
      source ${baseDir}/awk.sh 
      paste <(extract_fastqc -b ${fwd_pair}) <(extract_fastqc -o ${rvs_pair}) <(extract_fastqc -o ${fwd_unpair}) <(extract_fastqc -o ${rvs_unpair}) | body sort -k1 > qc_cleaned.txt
        
      """ 
  
     }
     
    }
}

/*
 * Create a QC summary statistics table 
 * Input: [samplename]_{1,2}P_Trimmed.txt , [samplename]_{1,2}U_Trimmed.txt
 * Output: "qc_pair.txt, qc_unpair.txt"
 */

process CREATE_QCTABLE {
  publishDir "$params.output/qc_table", mode: 'copy'
  input:
  path qc_raw
  path qc_cleaned
  path trim_sum

  output:
  path "*qc_table.txt" 

  script:
  if (params.SE){
  """
  #!/bin/bash 
  paste qc_raw.txt qc_cleaned.txt | \
  awk -F '\\t' '
  {
      \$18=sprintf("%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s",\$14,\$15,\$16,\$17,\$14,\$15,\$16,\$17)
  }
  NR==1{
      \$(NF+1)="Drop";
      } 
  NR>1{
      \$10= \$10 " (" \$10/\$2*100 ")";
      \$(NF+1)=\$2+\$6-\$10-\$19 " (" (\$2+\$6-\$10)/(\$2+\$6)*100 ")";
      } 
      {print}' OFS='\\t' > qc_table_tmp.txt
  echo -e "BEFORE_CLEANING\t\t\t\t\t\t\t\t\tAFTER_CLEANING" > header
  echo -e "FWD\t\t\t\t\tRVS\t\t\t\tPAIRED_FWD\t\t\t\tPAIRED_RVS\t\t\t\tUNPAIRED" >> header
  cat header qc_table_tmp.txt > qc_table.txt
  
  """

  }
  else if (params.mergeUnpaired){
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
      {print}' OFS='\\t' > qc_table_tmp.txt
  echo -e "BEFORE_CLEANING\t\t\t\t\t\t\t\t\tAFTER_CLEANING" > header
  echo -e "FWD\t\t\t\t\tRVS\t\t\t\tPAIRED_FWD\t\t\t\tPAIRED_RVS\t\t\t\tUNPAIRED" >> header
  cat header qc_table_tmp.txt > qc_table.txt
  
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
  echo -e "BEFORE_CLEANING\t\t\t\t\t\t\t\t\tAFTER_CLEANING" > header
  echo -e "FWD\t\t\t\t\tRVS\t\t\t\tPAIRED_FWD\t\t\t\tPAIRED_RVS\t\t\t\tUNPAIRED_FWD\t\t\t\tUNPAIRED_RVS" >> header
  cat header qc_table_tmp.txt > qc_table.txt
  """
  }
}

/*
 * Summarise FastQC results after cleaning using MultiQC
 * Input: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 * Output: "cleaned_data/ , cleaned.html"
 */

process MULTIQC_AFTER_CLEANING {
  publishDir "$params.output/multiqc", mode: 'copy'
  // errorStrategy 'ignore'

  input:
    path fastqc_after_all

  output:
    path "cleaned.html"
    path "cleaned_data"

  """
  multiqc ${fastqc_after_all} --interactive --filename cleaned 
  """
}

/* *** DEPRECATED process ***
 * Merge extracted QC results from before and after cleaning into one table
 * Input: qc_table_raw.txt, qc_pair.txt, qc_unpair.txt
 * Output: qc_table.txt
 */

process MERGE_QCTABLE {
  publishDir "$params.output/qc_table", mode: 'copy'

  input:
    path raw_qc_table
    path cleaned_qc_table_pair
    path cleaned_qc_table_unpair

  output:
    path "*qc_table.txt"
  
  """ 
  timestamp=\$(date '+%d%m%y_%H%M%S')
  paste ${raw_qc_table} ${cleaned_qc_table_pair} ${cleaned_qc_table_unpair} > \${timestamp}_qc_table.txt
  """
}

/* *** DEPRECATED process ***
 * Merge all samples' FastQC results into one table
 * Input: [samplename]_raw.txt
 * Output: qc_table_raw.txt
 */

process CREATE_QCTABLE_BEFORE_TRIM {
  publishDir "$params.output/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "qc_raw.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "${params.qc_header}") qc_table > qc_raw.txt
  """
}
