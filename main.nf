/*
 * Copyright (c) 2020
 * Pornchai Matangkasombut Center for Microbial Genomics, Department of Microbiology, Faculty of Science, Mahidol University
 * 
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
 * 
 * This Source Code Form is "Incompatible With Secondary Licenses", as
 * defined by the Mozilla Public License, v. 2.0.
 */



/* 
 * Enable DSL 2 syntax
 */
nextflow.enable.dsl = 2

/*
 * Define default parameter values
 */

// paths & inputs
baseDir         = "$HOME"
// Input fastq in glob pattern format 
params.input    = "$PWD/*[12].fastq.gz"
// Output directory
params.output   = "$PWD/output"

// Read type (true is PE)
params.SE = false 

// Merge unpaired F and R reads' stats (PE only)
// params.mergeUnpaired
params.mergeUnpaired = false

// cleaning options following Trimmomatic option format
// params.qc_options **
params.qc_options = "MINLEN:70"

// dev
/* params.qc_header = "Filename\tTotalSeq\tPoorQualSeq\tLength\t%GC\tavgSeqQual(min,max)" */
/* params.trim_header = "BothSurvied\tForwardOnlySurvived\tReverseOnlySurvied\tDroppedRead" */
time = new Date()

log.info """\
========================================================
DATE : $time
Path to sequence reads: $params.input
Single-end reads: $params.SE
Output directory : $params.output
Trimmomatic cleaning options: $params.qc_options
Unpaired reads' stats merged: $params.mergeUnpaired
========================================================
"""


/* 
 * Import modules
 * from ./modules.nf
 */
include { 
  FASTQC_BEFORE_CLEANING_SE;
  FASTQC_BEFORE_CLEANING_PE;
  EXTRACT_FASTQC_BEFORE_CLEANING;
  MULTIQC_BEFORE_CLEANING;
  CLEAN_SE;
  CLEAN_PE;
  EXTRACT_TRIM_LOG;
  CREATE_TRIM_SUMMARY_TABLE;
  FASTQC_AFTER_CLEANING_SE;
  FASTQC_AFTER_CLEANING_PE;
  EXTRACT_FASTQC_AFTER_CLEANING;
  CREATE_QCTABLE;
  MULTIQC_AFTER_CLEANING;
  MERGE_UNPAIRED_READS;
} from './modules.nf' 


/* 
 * Main pipeline 
 */
 
workflow {
  // input: single/paired-end reads
  // STEP 1: Quality check before cleaning
  if (params.SE) {
      reads = Channel.fromPath(params.input)
      FASTQC_BEFORE_CLEANING_SE(reads).set { FASTQC_BEFORE_CLEANING }
  } else {
      reads = Channel.fromFilePairs(params.input)
      FASTQC_BEFORE_CLEANING_PE(reads).set { FASTQC_BEFORE_CLEANING }
  }

  MULTIQC_BEFORE_CLEANING(FASTQC_BEFORE_CLEANING.collect())
  EXTRACT_FASTQC_BEFORE_CLEANING(FASTQC_BEFORE_CLEANING.flatten().collect())

  // STEP 2a: SE read cleaning and QCing
  if (params.SE){
      CLEAN_SE(reads, params.qc_options).set { CLEAN }
      
      // extract TRIM log 
      EXTRACT_TRIM_LOG(CLEAN_SE.out[1])
      CREATE_TRIM_SUMMARY_TABLE(EXTRACT_TRIM_LOG.out.collect())
      
      // QC cleaned reads 
      FASTQC_AFTER_CLEANING_SE(CLEAN_SE.out[0]).set { FASTQC_AFTER_CLEANING }

  } else {
  // STEP 2b: PE read cleaning and QCing
      CLEAN_PE(reads, params.qc_options).set { CLEAN }
      EXTRACT_TRIM_LOG(CLEAN_PE.out[2])
      CREATE_TRIM_SUMMARY_TABLE(EXTRACT_TRIM_LOG.out.collect())
      
      // QC cleaned reads
      // Merge unpaired F and R reads' stats (PE only)
      if(params.mergeUnpaired) {
          //MERGE_UNPAIRED_READS
        MERGE_UNPAIRED_READS(CLEAN_PE.out[1])
        FASTQC_AFTER_CLEANING_PE(
          CLEAN_PE.out[0],
          MERGE_UNPAIRED_READS.out
          ).set { FASTQC_AFTER_CLEANING }
      } else {
        FASTQC_AFTER_CLEANING_PE(
          CLEAN_PE.out[0],
          CLEAN_PE.out[1]
        ).set { FASTQC_AFTER_CLEANING }
      }
  }

  // STEP 3: Pooling FastQC results together using MultiQC 
  MULTIQC_AFTER_CLEANING(FASTQC_AFTER_CLEANING.collect())

  // STEP 4: Quality check after cleaning
  EXTRACT_FASTQC_AFTER_CLEANING(FASTQC_AFTER_CLEANING.flatten().collect())

  // STEP 5: Summarising the results
  CREATE_QCTABLE(
    EXTRACT_FASTQC_BEFORE_CLEANING.out,
    EXTRACT_FASTQC_AFTER_CLEANING.out,
    CREATE_TRIM_SUMMARY_TABLE.out
  )
  
}

  


