/*
 * Copyright (c) 2020
 * Center for Microbial Genomics, Mahidol University
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
 * Define the default parameters
 */

// paths & inputs
baseDir          = "$HOME"
// params.input     = "$baseDir/workspace/sumqc/Bp_fq/*_[1,2].[fq,fastq]*"
params.input     = "$baseDir/github/sumqc/test/input/test*_[1,2].fq.gz"
// params.results   = "$baseDir/workspace/sumqc/results"
params.results   = "$baseDir/github/sumqc/test/output"
// parameters for trimming
// params.qc_option = "SLIDINGWINDOW:4:30 MINLEN:70"
params.qc_option = "MINLEN:70"
params.mergeUnpair = false
params.SE = false 

params.qc_header = "Filename\tTotalSeq\tPoorQualSeq\tLength\t%GC\tavgSeqQual(min,max)"
params.trim_header = "BothSurvied\tForwardOnlySurvived\tReverseOnlySurvied\tDroppedRead"

log.info """\
========================================================
Raw reads: $params.input
Results directory : $params.results
Trimming option (for trimmomatic): $params.qc_option
========================================================
"""


/* 
 * Import modules 
 */
include { 
  FASTQC_BEFORE_TRIM_SE;
  FASTQC_BEFORE_TRIM_PE;
  EXTRACT_FASTQC_BEFORE_TRIM;
  /* CREATE_QCTABLE_BEFORE_TRIM; */
  MULTIQC_FASTQC_BEFORE_TRIM;
  TRIM_SE;
  TRIM_PE;
  EXTRACT_TRIM_LOG;
  CREATE_TRIM_SUMMARY_TABLE;
  FASTQC_AFTER_TRIM_SE;
  FASTQC_AFTER_TRIM_PE;
  EXTRACT_FASTQC_AFTER_TRIM;
  FINALISE_QCTABLE;
  MULTIQC_FASTQC_AFTER_TRIM;
  /* MERGE_QCTABLE; */
  MERGE_UNPAIRED_READ;
} from './modules.nf' 


/* 
 * Main pipeline steps
 */
 
workflow {
  // input: paired-end reads
    if (params.SE) {
        reads = Channel.fromPath(params.input)
        FASTQC_BEFORE_TRIM_SE(reads).set { FASTQC_BEFORE_TRIM }
    } else {
        reads = Channel.fromFilePairs(params.input)
        FASTQC_BEFORE_TRIM_PE(reads).set { FASTQC_BEFORE_TRIM }
    }

  // STEP 1: Quality check before cleaning
  EXTRACT_FASTQC_BEFORE_TRIM(FASTQC_BEFORE_TRIM.flatten().collect())

  // STEP 2: Cleaning
  if (params.SE){
      TRIM_SE(reads, params.qc_option).set { TRIM }
      EXTRACT_TRIM_LOG(TRIM_SE.out[1])
      CREATE_TRIM_SUMMARY_TABLE(EXTRACT_TRIM_LOG.out.collect())
      FASTQC_AFTER_TRIM_SE(TRIM_SE.out[0]).set { FASTQC_AFTER_TRIM }
  } else {
      TRIM_PE(reads, params.qc_option).set { TRIM }
      EXTRACT_TRIM_LOG(TRIM_PE.out[2])
      CREATE_TRIM_SUMMARY_TABLE(EXTRACT_TRIM_LOG.out.collect())
      
      // STEP 3: Quality check after cleaning
      if(params.mergeUnpair) {
        MERGE_UNPAIRED_READ(TRIM_PE.out[1])
        FASTQC_AFTER_TRIM_PE(
          TRIM_PE.out[0],
          MERGE_UNPAIRED_READ.out
          ).set { FASTQC_AFTER_TRIM }

      } else{
        FASTQC_AFTER_TRIM_PE(
          TRIM_PE.out[0],
          TRIM_PE.out[1]
        ).set { FASTQC_AFTER_TRIM }
      }
  }
  EXTRACT_FASTQC_AFTER_TRIM(
    FASTQC_AFTER_TRIM.flatten().collect()
  )
  FINALISE_QCTABLE(
    EXTRACT_FASTQC_BEFORE_TRIM.out,
    EXTRACT_FASTQC_AFTER_TRIM.out,
    CREATE_TRIM_SUMMARY_TABLE.out
  )
  MULTIQC_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.collect())
  
}

  


