/*
 * Process definitions
 */

/*
 * Step 1a: Read QC before trimming using fastqc
 * Input: Raw reads
 * Output: [samplename].fastqc.zip
 */

process FASTQC_BEFORE_TRIM {
  publishDir "$params.results/fastqc", mode: 'copy'
  tag "$id"

  input:
    tuple val(id), path(reads)

  output:
    path "*.zip"

  """
  fastqc $reads --noextract --quiet
  """
}

/*
 * Step 1b: Extract QC data of each sample 
 * Input: [samplename].fastqc.zip
 * Output: [samplename]_raw.txt
 */

process EXTRACT_FASTQC_BEFORE_TRIM {
  input: 
    path input

  output:
    path "*_raw.txt" 

  
  shell:
  '''
  unzip -p !{input} !{input.baseName}/fastqc_data.txt > !{input.baseName}

  awk -F '\\t' '/Filename|Total Sequences|Sequences flagged as poor quality|Sequence length|%GC/ {
        print $2
        }' !{input.baseName} > !{input.baseName}.txt \

  sed -n '/Per base sequence quality/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' ORS=" ">> !{input.baseName}.txt \

  sed -n '/Per sequence quality scores/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR==3 {
          print "(" $1 ", "
      } ORS=""; 
      {
          max=lastline;lastline=$1
      }; 
      END {
          print max ")\\n"
      }'  >> !{input.baseName}.txt \

  sed -n '/Per base N content/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' >> !{input.baseName}.txt \

  sed -n '/Total Deduplicated Percentage/,/END/p' !{input.baseName} |
  awk -F '\\t' 'NR==1 {print $2}' >> !{input.baseName}.txt \


  paste -s !{input.baseName}.txt > !{input.baseName}_raw.txt
  
  '''
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
  path "qc_table_raw.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "${params.qc_header}") qc_table > qc_table_raw.txt
  """


}
/*
 * Step 1d: Summarise fastqc results using multiqc
 * Input: [samplename].fastqc.zip
 * Output: "raw_data/ , raw.html"
 */
process MULTIQC_FASTQC_BEFORE_TRIM {
  publishDir "$params.results/multiqc_fastqc", mode: 'copy'
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
  errorStrategy 'ignore'

  input:
    tuple val(id), path(reads)
    path adapter
    val trimming_option

  output:
    tuple val(id), path("*.fq.gz")
    path "*_summary.txt"
    

  """
    TrimmomaticPE -trimlog ${id}_log.txt -summary ${id}_summary.txt\
    -basein ${reads[0]} \
    -baseout ${id}.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 ${trimming_option} 
  """
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
  path "trim_log.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "${params.trim_header}") qc_table > trim_log.txt

  """
}

/*
 * Step 3a: Read QC after trimming using fastqc
 * Input: [samplename]_{1,2}P.fq.gz, [samplename]_{1,2}U.fq.gz
 * Output: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 */
process FASTQC_AFTER_TRIM {
  tag "$id"

  input:
    tuple val(id), path(reads)

  output:
    path "*.zip"

  """
  fastqc $reads --noextract --quiet
  """
}

/*
 * Step 3b: Extract QC data from cleaned samples 
 * Input: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 * Output: [samplename]_{1,2}P_Trimmed.txt , [samplename]_{1,2}U_Trimmed.txt
 */

process EXTRACT_FASTQC_AFTER_TRIM {
  input: 
    path input

  output:
    path "*_fastqc_Trimmed.txt" 
  
  shell:
  '''
  unzip -p !{input} !{input.baseName}/fastqc_data.txt > !{input.baseName}

  awk -F '\\t' '/Filename|Total Sequences|Sequences flagged as poor quality|Sequence length|%GC/ {
        print $2
        }' !{input.baseName} > !{input.baseName}.txt \

  sed -n '/Per base sequence quality/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' ORS=" ">> !{input.baseName}.txt \

  sed -n '/Per sequence quality scores/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR==3 {
          print "(" $1 ", "
      } ORS=""; 
      {
          max=lastline;lastline=$1
      }; 
      END {
          print max ")\\n"
      }'  >> !{input.baseName}.txt \

  sed -n '/Per base N content/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' >> !{input.baseName}.txt \

  sed -n '/Total Deduplicated Percentage/,/END/p' !{input.baseName} |
  awk -F '\\t' 'NR==1 {print $2}' >> !{input.baseName}.txt \


  paste -s !{input.baseName}.txt > !{input.baseName}_Trimmed.txt
  
  '''
}

/*
 * Step 3c: Merge QC data of paired and unpaired samples into two different table
 * Input: [samplename]_{1,2}P_Trimmed.txt , [samplename]_{1,2}U_Trimmed.txt
 * Output: "qc_pair.txt, qc_unpair.txt"
 */

process CREATE_QCTABLE_AFTER_TRIM {
  publishDir "$params.results/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "*.txt" 

  """
  cat ${all_data} | sort > qc_table
  grep P.fq.gz qc_table > qc_pair
  grep U.fq.gz qc_table > qc_unpair
  cat <(echo "${params.qc_header}") qc_pair > qc_pair.txt
  cat <(echo "${params.qc_header}") qc_unpair  > qc_unpair.txt

  """
}

/*
 * Step 3d: Summarise fastqc results using multiqc
 * Input: [samplename]_{1,2}P.fastqc.zip , [samplename]_{1,2}U.fastqc.zip
 * Output: "cleaned_data/ , cleaned.html"
 */

process MULTIQC_FASTQC_AFTER_TRIM {
  publishDir "$params.results/multiqc_fastqc", mode: 'copy'
  errorStrategy 'ignore'

  input:
    path fastqc_after_all

  output:
    path "cleaned.html"
    path "cleaned_data"

  """
  multiqc ${fastqc_after_all} --interactive --filename cleaned
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
    path trim_qc_table

  output:
    path "qc_table.txt"
    
  shell:
  '''
  awk 'BEGIN {
        OFS="\\t"
        } 
    #print header
      NR==1 {
        print $0,$0,$0,"Dropped"
        } 
    #save total read of raw to array a
      FNR>1&&NR==FNR {
        raw[FNR]=$0;
        a[FNR]=$2;
        next;
        }
    #add ratio of survived pair read
      FNR>1 && FILENAME=="qc_pair.txt" {
        pairSeq[FNR]=$2; 
        $2=$2 " (" $2/a[FNR] ")";
        b[FNR]=$0;
        next;
        } 
    #add ratio of survived unpair read
      FNR>1 && FILENAME=="qc_unpair.txt" {
        $(NF+1)=(a[FNR]-pairSeq[FNR]-$2) " (" (a[FNR]-pairSeq[FNR]-$2)/a[FNR] ")"; #add number and ratio of dropped read
        $2=$2 " (" $2/a[FNR] ")"; 
        print raw[FNR],b[FNR],$0; #join tables
        }' FS='\\t' OFS='\\t' !{raw_qc_table} !{trim_qc_table} > qc_table.txt
  
  '''
}
