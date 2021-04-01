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
params.reads     = "$baseDir/workspace/sumqc/Bp_fq/*_{1,2}.fq.gz"
params.results   = "$baseDir/workspace/sumqc/results"
// parameters for trimming
params.trimming_option = "SLIDINGWINDOW:4:30 MINLEN:70"
params.mergeUnpair = false

params.qc_header = "Filename\tTotalSeq\tPoorQualSeq\tLength\t%GC\tavgSeqQual(min,max)"
params.trim_header = "BothSurvied\tForwardOnlySurvived\tReverseOnlySurvied\tDroppedRead"

log.info """\
snpplet v0.2
========================================================
Raw reads: $params.reads
Results directory : $params.results
Trimming option (for trimmomatic): $params.trimming_option
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
  CREATE_QCTABLE_AFTER_TRIM;
  MULTIQC_FASTQC_AFTER_TRIM;
  MERGE_QCTABLE;
  MERGE_UNPAIRED_READ;
} from './modules.nf' 


/* 
 * Main pipeline steps
 */
workflow {
  // input: paired-end reads
  read_pairs = Channel.fromFilePairs(params.reads)

  // STEP 1: Quality check before cleaning
  FASTQC_BEFORE_TRIM(read_pairs)
  EXTRACT_FASTQC_BEFORE_TRIM(FASTQC_BEFORE_TRIM.out.flatten())
  CREATE_QCTABLE_BEFORE_TRIM(EXTRACT_FASTQC_BEFORE_TRIM.out.collect())
  MULTIQC_FASTQC_BEFORE_TRIM(FASTQC_BEFORE_TRIM.out.collect())

  // STEP 2: Cleaning
  TRIM(
    read_pairs,
    params.trimming_option
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
    EXTRACT_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.flatten())
    CREATE_QCTABLE_AFTER_TRIM(EXTRACT_FASTQC_AFTER_TRIM.out.collect(),params.mergeUnpair)
    MULTIQC_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.collect())
    MERGE_QCTABLE(
    CREATE_QCTABLE_BEFORE_TRIM.out,
    CREATE_QCTABLE_AFTER_TRIM.out
    )
  }
  else{
    FASTQC_AFTER_TRIM(TRIM.out[0],TRIM.out[1])
    EXTRACT_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.flatten())
    CREATE_QCTABLE_AFTER_TRIM(EXTRACT_FASTQC_AFTER_TRIM.out.collect(),params.mergeUnpair)
    MULTIQC_FASTQC_AFTER_TRIM(FASTQC_AFTER_TRIM.out.collect())
    MERGE_QCTABLE(
    CREATE_QCTABLE_BEFORE_TRIM.out,
    CREATE_QCTABLE_AFTER_TRIM.out
    )
  }

  // STEP 4: Get all data from step 1 and 3 into single table  
  

}
