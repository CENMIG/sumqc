# sumqc

Goals:

- [x] Output the summary of reads quality before and after trimming into tab-delimited format

Simplified workflow:
1) Quality check raw reads [*fastqc*]
2) Trimming [*trimmomatic*]
3) Quality check again [*fastqc*]
4) Summarise statistic from 1 and 3 into one table [*bashscript*]

Pipelines are managed by *nextflow*. Most of the codes are reused and modified from snpplet.
