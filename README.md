# sumQC

  sumQC contains the [Nextflow](https://www.nextflow.io/) pipeline for quality checking and cleaning of both short reads and long reads sequencing data and also produce the summary of quality before and after cleaning into one table. The following software are use in the pipeline. 
 
|Pipeline|QC|Cleaning|
|:--:|:--:|:--:|
|Short reads|FastQC|Trimmomatic|
|Long reads|LongQC|Filtlong|

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

All related programs in pipeline are already included but their dependencies still need to install manually.
* linux OS
* java 8 or later
  * for [Nextflow](https://www.nextflow.io/), [FastQC](https://github.com/s-andrews/FastQC), and [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
* python3 
  * for LongQC, please check its dependencies [here](https://github.com/yfukasawa/LongQC)
* C++ complier (GCC 4.8) 
  * for [FiltLong](https://github.com/rrwick/Filtlong)


***

# Usage 

```
# For short reads
./sumqc SR -t <SE|PE> -i <path/to/shortread.fastq.gz> -o <path/for/output> -q <trimmomatic cleaning option>

# For long reads
./sumqc LR -i <path/to/longread.fastq.gz> -o <path/for/output> -q <FiltLong cleaning option>
```
If cleaning option is not specify the default is <SLIDINGWINDOW:4:30 MINLEN:70> for SR pipeline and <--min_length 10000 --keep_percent 90> for LR pipeline. 
Please see more details on cleaning option are available on their webpage ([Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [FiltLong](https://github.com/rrwick/Filtlong)).

Other available options are: 

`-h`      Print help
  
`-i`      Input directory

`-o`      Output directory

`-m`      Reserve for SR pipeline. this will merge unpair reads produced by Trimmomatic into single file.

`-q`      QC option (default:SR:SLIDINGWINDOW:4:30 MINLEN:70 ,LR:length 10000 --keep_percent 90)

`-r`      This give argument -resume to Nextflow which allow the pipeline to continue from where it stop in the previous run" 

