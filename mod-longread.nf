
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
         END{asort(a); printf "${raw_id.baseName} %s %.0f (%s-%s) %.2f %.2f\\n", num, n50, a[1], a[FNR], gc*100, sum/sumn}' FS=': ' ${qc_raw}/QC_vals_longQC_sampleqc.json FS='\\t' ${qc_raw}/longqc_sdust.txt > ${raw_id.baseName}_rawqc
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
    echo -e "filename\\tTotalSeq\\tN50(min,max)\\t%GC\\tavgSeqQual" | cat - $qc_raw > qc_raw.txt 
    
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
    input:
        path cleaned_id
        path qc_cleaned
    output:
        path "*_cleanedqc"
    shell:
    """
    awk 'NR==FNR&&/Num_of_reads/{gsub(",",""); num=\$2} 
         NR==FNR&&/N50/{n50=\$2} 
         NR==FNR&&/Mean_GC_content/{gc=\$2; next} 
         {a[NR]=\$3; \$7=\$3*\$5; ;sum+=\$7 ;sumn+=\$3} 
         END{asort(a); printf "${cleaned_id.baseName} %s %.0f (%s-%s) %.2f %.2f\\n", num, n50, a[1], a[FNR], gc*100, sum/sumn}' FS=': ' ${qc_cleaned}/QC_vals_longQC_sampleqc.json FS='\\t' ${qc_cleaned}/longqc_sdust.txt > ${cleaned_id.baseName}_cleanedqc
    """
}


//Concatenate to table
process CREATE_CLEAN_QC_TABLE {
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
    echo -e "filename\\tTotalSeq\\tN50(min,max)\\t%GC\\tavgSeqQual" | cat - $qc_cleaned > qc_cleaned.txt 
    paste $qc_raw qc_cleaned.txt > \${timestamp}qc_table.txt    
    """
}




