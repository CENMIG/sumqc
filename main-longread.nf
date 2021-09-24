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
 * Define the default parameter values
 */

// paths & inputs
baseDir          = "$HOME"
// Input fastq in glob pattern format 
params.input     = "$PWD/*.fastq.gz"
// Output directory
params.output    = "$PWD/output"

// Sequencing platform
params.platform  = "pb-rs2"

// Cleaning options following Trimmomatic option format 
params.qc_options = "--min_length 10000 --keep_percent 90"

time = new Date()

log.info """\
========================================================
DATE: $time
Path to sequence reads: $params.input
Output directory: $params.output
Sequencing platform: $params.platform
Filtlong cleaning options: $params.qc_options
========================================================
"""


/* 
 * Import modules 
 * from ./mod-longread.nf
 */

include { 
    LONGQC_BEFORE_CLEANING
    EXTRACT_LONGQC_BEFORE_CLEANING
    CREATE_QCTABLE_BEFORE_CLEANING
    CLEAN
    LONGQC_AFTER_CLEANING
    EXTRACT_LONGQC_AFTER_CLEANING
    CREATE_QCTABLE
} from './mod-longread.nf'


/* 
 * Main pipeline 
 */

workflow {
    // input: long reads fastq file
    inputChannel=Channel.fromPath(params.input)
    
    // Step 1: Input fastq and quality check before cleaning 
    LONGQC_BEFORE_CLEANING(inputChannel,params.platform)
    EXTRACT_LONGQC_BEFORE_CLEANING(inputChannel,LONGQC_BEFORE_CLEANING.out)
    CREATE_QCTABLE_BEFORE_CLEANING(EXTRACT_LONGQC_BEFORE_CLEANING.out.collect())

    // Step 2: Read cleaning 
    CLEAN(inputChannel,params.qc_options)

    // STEP 3: Quality check after cleaning
    LONGQC_AFTER_CLEANING(CLEAN.out,params.platform)
    EXTRACT_LONGQC_AFTER_CLEANING(CLEANING.out,LONGQC_AFTER_CLEANING.out)

    // STEP 4: Create QC table
    CREATE_QCTABLE(CREATE_QCTABLE_BEFORE_CLEANING.out,EXTRACT_LONGQC_AFTER_CLEANING.out.collect())
}
