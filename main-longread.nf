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
params.input     = "$baseDir/workspace/longread/PacBio_RS_II/*.fq.gz"
params.results   = "$baseDir/workspace/results"
params.longQC    = "$baseDir/programs/LongQC/longQC.py"
params.filtLong  = "$baseDir/programs/Filtlong/bin/filtlong"
// parameters for trimming
params.platform = "pb-rs2"
params.qc_option = "--min_length 10000 --keep_percent 90"
params.cleaning_option = "SLIDINGWINDOW:4:30 MINLEN:70"

//params.qc_header = "Filename\tTotalSeq\tLength\t%GC\tavgSeqQual(min,max)"
//params.trim_header = "BothSurvied\tForwardOnlySurvived\tReverseOnlySurvied\tDroppedRead"

log.info """\
snpplet v0.2
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
    QC_RAW_DATA
    EXTRACT_RAW_QC
    CREATE_RAW_QC_TABLE
    CLEANING
    QC_CLEAN_DATA
    EXTRACT_CLEAN_QC
    CREATE_CLEAN_QC_TABLE
} from './mod-longread.nf'


/* 
 * Main pipeline steps
 */
workflow {
    inputChannel=Channel.fromPath(params.input)
    QC_RAW_DATA(inputChannel,params.platform)
    EXTRACT_RAW_QC(inputChannel,QC_RAW_DATA.out)
    CREATE_RAW_QC_TABLE(EXTRACT_RAW_QC.out.collect())
    CLEANING(inputChannel,params.qc_option)
    QC_CLEAN_DATA(CLEANING.out,params.platform)
    EXTRACT_CLEAN_QC(CLEANING.out,QC_CLEAN_DATA.out)
    CREATE_CLEAN_QC_TABLE(EXTRACT_CLEAN_QC.out.collect())
}
