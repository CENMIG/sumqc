# sumQC

sumQC contains the pipeline for quality checking and cleaning of both short reads and long reads sequencing data and also produce the summary of quality before and after cleaning into one table. The following software are use in the pipeline.

||QC|Cleangin|
|--|--|--|
|SR|FastQC|Trimmomatic|
|LR|LongQC|Filtlong|

Simplified workflow:
1) Quality check raw reads. [*fastqc*/*LongQC*]
2) Cleaning. [*trimmomatic*/*FiltLong*]
3) Quality check again. [*fastqc*/*LongQC*]
4) Summarise quality statistic into one table.
***

# Installation 

```
git clone https://github.com/thavinb/sumqc.git
./sumqc -h
```
The program come with wrapper name 'sumqc'. Consider put it in your PATH for convenient usage.

***

# Requirement

All related programs use in pipeline also come inside the prog directory. Both program for QC and cleaning for short reads are written in java. While program for the long reads are written in python. 
***

# Usage 

```
# For short reads
./sumqc SR -i <path/to/shortread.fastq.gz> -o <path/for/output> -q <trimmomatic cleaning option>
# For long reads
./sumqc LR -i <path/to/longread.fastq.gz> -o <path/for/output> -q <FiltLong cleaning option>
```




Goals:

- [x] Output the summary of reads quality before and after trimming into tab-delimited format


Pipelines are managed by *nextflow*. Most of the codes are reused and modified from snpplet.
