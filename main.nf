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
 * 'snpplet' - A Nextflow pipeline for variant calling from NGS data
 * 
 * Yuttapong Thawornwattana
 * Bharkbhoom Jaemsai
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
params.input     = "$baseDir/workspace/sumqc/Bp_fq/TESTA_[1,2].[fq,fastq]*"
// params.results   = "$baseDir/workspace/sumqc/results"
params.results   = "$baseDir/workspace/sumqc/dev"
// parameters for trimming
// params.qc_option = "SLIDINGWINDOW:4:30 MINLEN:70"
params.qc_option = "MINLEN:70"
params.mergeUnpair = false

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
  FASTQC_BEFORE_TRIM;
  EXTRACT_FASTQC_BEFORE_TRIM;
  CREATE_QCTABLE_BEFORE_TRIM;
  MULTIQC_FASTQC_BEFORE_TRIM;
  TRIM;
  EXTRACT_TRIM_LOG;
  CREATE_TRIM_SUMMARY_TABLE;
  FASTQC_AFTER_TRIM;
  EXTRACT_FASTQC_AFTER_TRIM;
  FINALISE_QCTABLE;
  MULTIQC_FASTQC_AFTER_TRIM;
  MERGE_QCTABLE;
  MERGE_UNPAIRED_READ;
} from './modules.nf' 


/* 
 * Main pipeline steps
 */
workflow {
  // input: paired-end reads
  read_pairs = Channel.fromFilePairs(params.input)

  // STEP 1: Quality check before cleaning
  FASTQC_BEFORE_TRIM(read_pairs)
  FWD_RAW=FASTQC_BEFORE_TRIM.out.flatten().filter{ it =~ /.*_1_fastqc.zip$/ }
  RVS_RAW=FASTQC_BEFORE_TRIM.out.flatten().filter{ it =~ /.*_2_fastqc.zip$/ }
  EXTRACT_FASTQC_BEFORE_TRIM(FWD_RAW.collect(),RVS_RAW.collect())
  // CREATE_QCTABLE_BEFORE_TRIM(EXTRACT_FASTQC_BEFORE_TRIM.out.collect())
  // MULTIQC_FASTQC_BEFORE_TRIM(FASTQC_BEFORE_TRIM.out.collect())

  // STEP 2: Cleaning
  TRIM(
    read_pairs,
    params.qc_option
  )
  EXTRACT_TRIM_LOG(TRIM.out[2])
  CREATE_TRIM_SUMMARY_TABLE(EXTRACT_TRIM_LOG.out.collect())
  
  // STEP 3: Quality check after cleaning
  if(params.mergeUnpair) {
    MERGE_UNPAIRED_READ(TRIM.out[1])
    FASTQC_AFTER_TRIM(
      TRIM.out[0],
      MERGE_UNPAIRED_READ.out
      )
    PFWD_CLEANED=FASTQC_AFTER_TRIM.out.flatten().filter{ it =~ /.*_1P_fastqc.zip$/ }
    PRVS_CLEANED=FASTQC_AFTER_TRIM.out.flatten().filter{ it =~ /.*_2P_fastqc.zip$/ }
    UFWD_CLEANED=FASTQC_AFTER_TRIM.out.flatten().filter{ it =~ /.*_U_fastqc.zip$/ }
    EXTRACT_FASTQC_AFTER_TRIM(
        PFWD_CLEANED.collect(),
        PRVS_CLEANED.collect(),
        UFWD_CLEANED.collect(),
        URVS_CLEANED.collect())
    FINALISE_QCTABLE(
        EXTRACT_FASTQC_BEFORE_TRIM.out,
        EXTRACT_FASTQC_AFTER_TRIM.out,
        CREATE_TRIM_SUMMARY_TABLE.out,
        params.mergeUnpair)
    MULTIQC_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.collect())
    //MERGE_QCTABLE(
    //CREATE_QCTABLE_BEFORE_TRIM.out,
    //CREATE_QCTABLE_AFTER_TRIM.out[0],
    //CREATE_QCTABLE_AFTER_TRIM.out[1],
   // )
  }
  else{
    FASTQC_AFTER_TRIM(TRIM.out[0],TRIM.out[1])
    PFWD_CLEANED=FASTQC_AFTER_TRIM.out.flatten().filter{ it =~ /.*_1P_fastqc.zip$/ }
    PRVS_CLEANED=FASTQC_AFTER_TRIM.out.flatten().filter{ it =~ /.*_2P_fastqc.zip$/ }
    UFWD_CLEANED=FASTQC_AFTER_TRIM.out.flatten().filter{ it =~ /.*_1U_fastqc.zip$/ }
    URVS_CLEANED=FASTQC_AFTER_TRIM.out.flatten().filter{ it =~ /.*_2U_fastqc.zip$/ }
    EXTRACT_FASTQC_AFTER_TRIM(
        PFWD_CLEANED.collect(),
        PRVS_CLEANED.collect(),
        UFWD_CLEANED.collect(),
        URVS_CLEANED.collect())
    FINALISE_QCTABLE(
        EXTRACT_FASTQC_BEFORE_TRIM.out,
        EXTRACT_FASTQC_AFTER_TRIM.out,
        CREATE_TRIM_SUMMARY_TABLE.out,
        params.mergeUnpair)
    MULTIQC_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.collect())
    //MERGE_QCTABLE(
    //CREATE_QCTABLE_BEFORE_TRIM.out,
    //CREATE_QCTABLE_AFTER_TRIM.out
    //)
  }

  // STEP 4: Get all data from step 1 and 3 into single table  
  

}
