# sumQC

  sumQC is the pipeline use for checking quality, cleaning, and making the table of QC statistic before and after cleaning. The pipeline is built using [Nextflow](https://www.nextflow.io/), a workflow tool to run tasks across multiple compute infrasturcture.
 
  sumQC contains two workflows for short read and long read which use two different sets of program for quality checking and cleaning. Workflow must be specified as a first argument before adding other option (SR for short read and LR for long read). The programs used in the workflows are: 
  
|Workflow|QC|Cleaning|
|:--:|:--:|:--:|
|Short read|[FastQC](https://github.com/s-andrews/FastQC)|[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)|
|Long read|[LongQC](https://github.com/yfukasawa/LongQC)|[FiltLong](https://github.com/rrwick/Filtlong)|

Cleaning options are use according to the program used in the workflow. Please refers to the  program webpage for the full list of options.

***

# Pipeline summary 

1) Quality check before cleaning. `fastqc` `multiqc` `LongQC`
2) Extract quality statistic before cleaning `bash script`
3) Cleaning fastq file. `Trimmomatic` `FiltLong`
4) Quality check before cleaning. `fastqc` `multiqc` `LongQC`
5) Extract quality statistic before cleaning `bash script`
6) Summarise quality statistic into one table.


***

# Requirement

Most of the programs in pipeline are alredy included.

1. Java 8 or late 
2. [Nextflow](https://www.nextflow.io/) (>=21.04.0)
3. [Multiqc](https://github.com/ewels/MultiQC.git)
```
# Install through pip
pip install multiqc

# Install through conda
conda install -c bioconda multiqc
```
4. Python library dependecy for LongQC (please check [LongQC github page](https://github.com/yfukasawa/LongQC) for more details)
```
# Install through pip
pip install numpy scipy matplotlib scikit-learn pandas jinja2 h5py pysam edlib python-edlib

# Install through conda
conda install h5py
conda install -c bioconda pysam
conda install -c bioconda edlib
conda install -c bioconda python-edlib

```

***

# Installation 

```
git clone https://github.com/thavinb/sumqc.git
./sumqc -h
```
The program come with wrapper name 'sumqc'. Consider put it in your PATH for convenient usage.

***

# Usage 

   sumQC takes a glob pattern of fastq file as input. If no pattern are specified the program will use default glob pattern according to the workflow and the type of read user specified (i.e. *.fastq.gz or *[12].fastq.gz).
    
```
# For short reads

./sumqc SR -t <SE|PE> -i <path/to/shortread.fastq.gz> -o <path/for/output> -q <trimmomatic cleaning option>

./sumqc SR -t SE -i 'path/to/shortreads/*.fastq.gz' -o path/to/output/ -q "MINLEN:70 AVGQUAL:30" 
or 
./sumqc SR -t PE -i 'path/to/longreads/*[12].fastq.gz' -o path/to/output/ -q "MINLEN:70 AVGQUAL:30" 

# For long reads

./sumqc LR -p <pb-rs2|pb-sequel|ont-ligation|ont-rapid|ont-1dsq> -i <path/to/longread.fastq.gz> -o <path/for/output> -q <FiltLong cleaning option>

./sumqc LR -p ont-rapid -i path/to/longread/*.fastq.gz -o path/to/output/ -q --min_length 10000 --keep_percent 90
```


All available options are: 
```
-h      Print this help

-i      Glob pattern of input fastq file
        (i.e. *[12].fastq.gz)

-o      Output directory

-m      For paired-end read. Merge unpaired 
        forward and reverse reads into one file 
        after cleaning process.

-p      Sequencing platform for long read.
        (https://github.com/yfukasawa/LongQC)

        Platforms           Options
        ---------           -------
        RS-II               pb-rs2
        Sequel              pb-sequel
        ONT(1D ligation)    ont-ligation
        ONT(rapid)          ont-rapid

-q      Cleaning options 

        default for short read workflow:
        SLIDINGWINDOW:4:30 MINLEN:70

        default for long read workflow:
        min_length 10000 --keep_percent 90)

-r      Resume from cache in work directoty.
        This is the same as:

        nextflow run main.nf -resume.

        (https://www.nextflow.io/docs/latest/cli.html?highlight=resume#run)
```
***

# Output 

sumQC will output all the results of programs used in the pipeline and the QC statistic table of fastq before and after cleaning. 
The format of the table shown below: 

|RAW|||||||||CLEANED|||||||||||||||||
|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|
|FWD|||||RVS||||PAIRED_FWD||||PAIRED_RVS||||UNPAIRED_FWD||||UNPAIRED_RVS|||||
|FileName|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|TotalSeq|Length|GC|AvgQScore (min,max)|Drop|
|test_1-1T|10000|150|65.000|33.828 (27,37)|10000|150|66.000|31.439 (26,37)|2636 (26.36)|70-150|62.000|35.872 (34,37)|2636 (26.36)|70-150|61.000|35.483 (33,37)|3796 (37.96)|70-150|65.000|35.606 (34,37)|819 (8.19)|70-150|61.000|35.258 (33,37)|10113 (50.565)|
|test_3-1T|10000|150|66.000|33.917 (27,37)|10000|150|66.000|31.542 (26,37)|2900 (29)|70-150|63.000|35.941 (34,37)|2900 (29)|70-150|62.000|35.490 (33,37)|3748 (37.48)|70-150|66.000|35.676 (34,37)|784 (7.84)|70-150|62.000|35.310 (33,37)|9668 (48.34)|

***

# Acknowledgements 

An extensive list of references for the tools used by the pipeline can be found in the CITATIONS.md file.


