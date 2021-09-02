# sumQC

  sumQC contains the [Nextflow](https://www.nextflow.io/) pipeline for quality checking and cleaning of both short reads and long reads sequencing data and also produce the summary of quality before and after cleaning into one table. The following software are use in the pipeline. 
 
|Pipeline|QC|Cleaning|
|:--:|:--:|:--:|
|Short reads|FastQC|Trimmomatic|
|Long reads|LongQC|Filtlong|

Simplified workflow:
1) Quality check raw reads. [*fastqc*/*multiqc*/*LongQC*]
2) Cleaning. [*trimmomatic*/*FiltLong*]
3) Quality check again. [*fastqc*/*multiqc*/*LongQC*]
4) Summarise quality statistic into one table.


***

# Requirement

All related programs in pipeline are already included but their dependencies still need to install manually.
* linux OS
* java 8 or later
  * for [Nextflow](https://www.nextflow.io/), [FastQC](https://github.com/s-andrews/FastQC), and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* python3 
  * for LongQC, please check its dependencies [here](https://github.com/yfukasawa/LongQC)
* C++ complier (GCC 4.8) 
  * for [FiltLong](https://github.com/rrwick/Filtlong)

***

# Installation 

```
git clone https://github.com/thavinb/sumqc.git
./sumqc -h
```
The program come with wrapper name 'sumqc'. Consider put it in your PATH for convenient usage.

***


# Usage 

```
# For short reads
./sumqc SR -t <SE|PE> -i <path/to/shortread.fastq.gz> -o <path/for/output> -q <trimmomatic cleaning option>
./sumqc SR -t SE -i path/to/shortreads/*.fq.gz -o path/to/output/ -q "MINLEN:70 AVGQUAL:30" 
or 
./sumqc SR -t PE -i path/to/longreads/*[12].fq.gz -o path/to/output/ -q "MINLEN:70 AVGQUAL:30" 

# For long reads
./sumqc LR -p <pb-rs2|pb-sequel|ont-ligation|ont-rapid|ont-1dsq> -i <path/to/longread.fastq.gz> -o <path/for/output> -q <FiltLong cleaning option>
./sumqc LR -p ont-rapid -i path/to/longread/*.fq.gz -o path/to/output/ -q --min_length 10000 --keep_percent 90
```
If cleaning option is not specify the default for short reads and long reads pipeline are:
```
For short reads:
SLIDINGWINDOW:4:30 MINLEN:70
For Long reads:
--min_length 10000 --keep_percent 90
```
Please to the programs's webpage for more details on the options. ([Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [FiltLong](https://github.com/rrwick/Filtlong)).

Other available options are: 

`-h`      Print help
  
`-i`      Input directory

`-o`      Output directory

`-p`      Only use in LR pipeline to specify the sequencing platform. 
          The options are `<pb-rs2|pb-sequel|ont-ligation|ont-rapid|ont-1dsq>`.

`-m`      Only use for paired-end data in SR pipeline. this will merge unpair reads produced by Trimmomatic into single file.

`-q`      QC option (default:SR:SLIDINGWINDOW:4:30 MINLEN:70 ,LR:min_length 10000 --keep_percent 90)

`-r`      Allow the pipeline to continue from  it stop in the previous run 

***

# Results 

Others than the cleaned fastq files. sumQC's pipeline also produce simple qc table in the following format:


|RAW|||||||||CLEANED|||||||||||||||||
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|FWD|||||RVS||||PAIRED_FWD||||PAIRED_RVS||||UNPAIRED_FWD||||UNPAIRED_RVS|||||
|FileName|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|Drop|
|test_1-1T|10000|150|65.000|33.828 (27,37)|10000|150|66.000|31.439 (26,37)|2636 (26.36)|70-150|62.000|35.872 (34,37)|2636 (26.36)|70-150|61.000|35.483 (33,37)|3796 (37.96)|70-150|65.000|35.606 (34,37)|819 (8.19)|70-150|61.000|35.258 (33,37)|10113 (50.565)|
|test_3-1T|10000|150|66.000|33.917 (27,37)|10000|150|66.000|31.542 (26,37)|2900 (29)|70-150|63.000|35.941 (34,37)|2900 (29)|70-150|62.000|35.490 (33,37)|3748 (37.48)|70-150|66.000|35.676 (34,37)|784 (7.84)|70-150|62.000|35.310 (33,37)|9668 (48.34)|


