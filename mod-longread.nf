
//QC the raw data 
process QC_RAW_DATA {
    publishDir "$params.results/longQC", mode: 'copy'
    input:
        path raw_fastq
        val platform
    output:
	  path qc_raw 

// path to longQC.py consider to put it in main folder?
// add option for -x argument (pb-rs2 pb-sequel, ont-ligation, ont-rapid, ont-1dsq)
    shell:
   ''' 
    /usr/bin/env python3 !{baseDir}/prog/LongQC/longQC.py sampleqc -x !{platform} -o qc_raw !{raw_fastq} 
   '''
}


// Extract the QC raw data
process EXTRACT_RAW_QC {
    stageInMode "symlink"
    input:
        path raw_id
        path qc_raw
    output:
        path "*_rawqc"

    shell:
    """
    awk 'NR==FNR&&/Num_of_reads/{gsub(",",""); num=\$2} 
         NR==FNR&&/N50/{n50=\$2} 
         NR==FNR&&/Mean_GC_content/{gc=\$2; next} 
         {a[NR]=\$3; \$7=\$3*\$5; ;sum+=\$7 ;sumn+=\$3} 
         END{asort(a); printf "${raw_id.SimpleName}\\t%s\\t%s-%s\\t%.2f\\t%.2f\\n", num,  a[1], a[FNR], gc*100, sum/sumn}' FS=': ' ${qc_raw}/QC_vals_longQC_sampleqc.json FS='\\t' ${qc_raw}/longqc_sdust.txt > ${raw_id.SimpleName}_rawqc
    """
}

//Concatenate to table
process CREATE_RAW_QC_TABLE {
    publishDir "$params.results/qc_table", mode: 'copy'
    input:
        path qc_raw
    output: 
        path "qc_raw.txt"  

    shell:
    """
    echo -e "FileName\\tTotalSeq\\tLength\\t%GC\\tAvgQScore" | cat - $qc_raw > qc_raw.txt 
    
    """
}
// Cleaning data
process CLEANING {
    publishDir "$params.results/cleaned_fastq", mode: 'copy'
    input:
        path raw_fastq
        val qc_option
    output:
        path "*.fq.gz"

    shell:
    ''' 
    !{baseDir}/prog/Filtlong/bin/filtlong !{qc_option} !{raw_fastq} | gzip > !{raw_fastq.baseName}_cleaned.fq.gz
    '''
}


// QC cleaned data
process QC_CLEAN_DATA {
    publishDir "$params.results/longQC", mode: 'copy'
    input:
        path cleaned_fastq
        val platform
    output:
        path qc_cleaned

    shell:
    ''' 
    /usr/bin/env python3 !{baseDir}/prog/LongQC/longQC.py sampleqc -x !{platform} -o qc_cleaned !{cleaned_fastq} 

    '''
}

// Extract cleaned QC
process EXTRACT_CLEAN_QC {
    stageInMode "symlink"
    input:
        path cleaned_id
        path qc_cleaned
    output:
        path "*_cleanedqc"
    shell:
    """
    #/bin/bash
    name="${cleaned_id.SimpleName}"
    simplename=\${name%%_cleaned*}
    awk -v simplename=\${simplename} 'NR==FNR&&/Num_of_reads/{gsub(",",""); num=\$2} 
         NR==FNR&&/N50/{n50=\$2} 
         NR==FNR&&/Mean_GC_content/{gc=\$2; next} 
         {a[NR]=\$3; \$7=\$3*\$5; ;sum+=\$7 ;sumn+=\$3} 
         END{asort(a); printf "%s\\t%s\\t%s-%s\\t%.2f\\t%.2f\\n",simplename, num, a[1], a[FNR], gc*100, sum/sumn}' FS=': '  ${qc_cleaned}/QC_vals_longQC_sampleqc.json FS='\\t'  ${qc_cleaned}/longqc_sdust.txt > ${cleaned_id.SimpleName}_cleanedqc
    """
}


//Concatenate to table
process CREATE_QC_TABLE {
    publishDir "$params.results/qc_table", mode: 'copy'
    input:
        path qc_raw
        path qc_cleaned
    output: 
        path "qc_cleaned.txt"
        path "*qc_table.txt"

    shell:
    """
    timestamp=\$(date '+%d%m%y_%H%M%S')
    echo -e "FileName\\tTotalSeq\\tLength\\t%GC\\tAvgQScore" | cat - $qc_cleaned > qc_cleaned.txt 
    paste $qc_raw <(awk -F '\\t' 'NR==1{print \$2,\$3,\$4,\$5} NR>=2{for(i=1;i<=NF;i++) gsub(".*","NA",\$i);} NR>1{print \$2,\$3,\$4,\$5}' OFS='\\t' $qc_raw) > qc_raw.tmp
    paste qc_cleaned.txt <(awk -F '\\t' 'NR==1{print \$2,\$3,\$4,\$5,\$2,\$3,\$4,\$5,\$2,\$3,\$4,\$5} NR>=2{for(i=1;i<=NF;i++) gsub(".*","NA",\$i);} NR>1{print \$2,\$3,\$4,\$5,\$2,\$3,\$4,\$5,\$2,\$3,\$4,\$5}' OFS='\\t' qc_cleaned.txt) > qc_cleaned.tmp
    join --header -t \$'\\t' qc_raw.tmp qc_cleaned.tmp |  
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
  echo -e "FWD\t\t\t\t\tRVS\t\t\t\tPAIRED_FWD\t\t\t\tPAIRED_RVS\t\t\t\tUNPAIRED_FWD\t\t\t\tUNPAIRED_RVS" >> header
  cat header qc_table_tmp.txt > qc_table.txt    
    """
}




