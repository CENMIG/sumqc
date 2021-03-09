/*
 * Process definitions
 */


/*
 * Step 0a: Create a FASTA ref_genome index (.fai) with samtools for GATK
 */
process PREPARE_GENOME_SAMTOOLS { 
  tag "$ref_genome.baseName"
 
  input: 
    path ref_genome
 
  output: 
    path "${ref_genome}.fai"
  
  script:
  """
  samtools faidx $ref_genome
  """
}


/*
 * Step 0b: Create a FASTA ref_genome index with bwa for bwa
 */
process PREPARE_GENOME_BWA { 
  tag "$ref_genome.baseName"
 
  input: 
    path ref_genome
 
  output: 
    path "${ref_genome}.*"
  
  script:
  """
  bwa index -a is $ref_genome
  """
}


/*
 * Step 0c: Create a FASTA ref_genome sequence dictionary with Picard for GATK
 */
process PREPARE_GENOME_PICARD {
  tag "$ref_genome.baseName"

  input:
    path ref_genome

  output:
    path "${ref_genome.baseName}.dict"

  script:
  """
  gatk CreateSequenceDictionary -R $ref_genome -O ${ref_genome.baseName}.dict
  """
}


/*
 * Step 1a: Read QC before trimming using fastqc
 */
process FASTQC_BEFORE_TRIM {
  publishDir "$params.results/fastqc", mode: 'copy'
  tag "$id"

  input:
    tuple val(id), path(reads)

  output:
    path "*.zip"

  """
  fastqc $reads --noextract --quiet
  """
}

process EXTRACT_FASTQC_BEFORE_TRIM {
  input: 
    path input

  output:
    path "*_Tr.txt" 

  
  shell:
  '''
  unzip -p !{input} !{input.baseName}/fastqc_data.txt > !{input.baseName}

  awk -F '\\t' '/Filename|Total Sequences|Sequences flagged as poor quality|Sequence length|%GC/ {
        print $2
        }' !{input.baseName} > !{input.baseName}.txt \

  sed -n '/Per base sequence quality/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' ORS=" ">> !{input.baseName}.txt \

  sed -n '/Per sequence quality scores/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR==3 {
          print "(" $1 ", "
      } ORS=""; 
      {
          max=lastline;lastline=$1
      }; 
      END {
          print max ")\\n"
      }'  >> !{input.baseName}.txt \

  sed -n '/Per base N content/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' >> !{input.baseName}.txt \

  sed -n '/Total Deduplicated Percentage/,/END/p' !{input.baseName} |
  awk -F '\\t' 'NR==1 {print $2}' >> !{input.baseName}.txt \


  paste -s !{input.baseName}.txt > !{input.baseName}_Tr.txt
  
  '''
}

process CREATE_QCTABLE_BEFORE_TRIM {
  publishDir "$params.results/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "qc_table_raw.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "${params.qc_header}") qc_table > qc_table_raw.txt
  """


}
/*
 * Step 1b: Summarise fastqc results using multiqc
 */
process MULTIQC_FASTQC_BEFORE_TRIM {
  publishDir "$params.results/multiqc_fastqc", mode: 'copy'
  errorStrategy 'ignore'

  input:
    path fastqc_before_all

  output:
    path "before.html"
    path "before_data"

  """
  multiqc ${fastqc_before_all} --interactive --filename before
  """
}


/*
 * Step 1c. Trim adapters and low quality reads
 */
process TRIM {
  publishDir "$params.results/trim_results", mode: 'copy'
  tag "$id"
  errorStrategy 'ignore'

  input:
    tuple val(id), path(reads)
    path adapter
    val trimming_option

  output:
    tuple val(id), path("*.fq.gz")
    path "*_summary.txt"
    

  """
    TrimmomaticPE -trimlog ${id}_log.txt -summary ${id}_summary.txt\
    -basein ${reads[0]} \
    -baseout ${id}.fq.gz \
    ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 \
    ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:30 MINLEN:70
  """
}

process EXTRACT_TRIM_LOG {
  input:
  path summarylog

  output:
  path "*_ext.txt"

  """
  cat ${summarylog} | awk -F ':' '{print \$2}' | paste -s | 
  awk  '{print \$3,\$5,\$7,\$9}' OFS='\\t' > ${summarylog.baseName}_ext.txt
  """
}
process CREATE_TRIM_SUMMARY_TABLE {
  publishDir "$params.results/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "trim_log.txt"

  """
  cat  ${all_data} | sort > qc_table
  cat <(echo "${params.trim_header}") qc_table > trim_log.txt

  """
}
/*
 * Step 1d: Read QC after trimming using fastqc
 */
process FASTQC_AFTER_TRIM {
  tag "$id"

  input:
    tuple val(id), path(reads)

  output:
    path "*.zip"

  """
  fastqc $reads --noextract --quiet
  """
}

process EXTRACT_FASTQC_AFTER_TRIM {
  input: 
    path input

  output:
    path "*_fastqc_Trimmed.txt" 
  
  shell:
  '''
  unzip -p !{input} !{input.baseName}/fastqc_data.txt > !{input.baseName}

  awk -F '\\t' '/Filename|Total Sequences|Sequences flagged as poor quality|Sequence length|%GC/ {
        print $2
        }' !{input.baseName} > !{input.baseName}.txt \

  sed -n '/Per base sequence quality/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' ORS=" ">> !{input.baseName}.txt \

  sed -n '/Per sequence quality scores/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR==3 {
          print "(" $1 ", "
      } ORS=""; 
      {
          max=lastline;lastline=$1
      }; 
      END {
          print max ")\\n"
      }'  >> !{input.baseName}.txt \

  sed -n '/Per base N content/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' >> !{input.baseName}.txt \

  sed -n '/Total Deduplicated Percentage/,/END/p' !{input.baseName} |
  awk -F '\\t' 'NR==1 {print $2}' >> !{input.baseName}.txt \


  paste -s !{input.baseName}.txt > !{input.baseName}_Trimmed.txt
  
  '''
}

process CREATE_QCTABLE_AFTER_TRIM {
  publishDir "$params.results/qc_table", mode: 'copy'
  input:
  path all_data

  output:
  path "*.txt" 

  """
  cat ${all_data} | sort > qc_table
  grep P.fq.gz qc_table > qc_pair
  grep U.fq.gz qc_table > qc_unpair
  cat <(echo "${params.qc_header}") qc_pair > qc_pair.txt
  cat <(echo "${params.qc_header}") qc_unpair  > qc_unpair.txt

  """
}
/*
 * Step 1e: Summarise fastqc results using multiqc
 */
process MULTIQC_FASTQC_AFTER_TRIM {
  publishDir "$params.results/multiqc_fastqc", mode: 'copy'
  errorStrategy 'ignore'

  input:
    path fastqc_after_all

  output:
    path "after.html"
    path "after_data"

  """
  multiqc ${fastqc_after_all} --interactive --filename after
  """
}

process MERGE_QCTABLE {
  publishDir "$params.results/qc_table", mode: 'copy'

  input:
    path raw_qc_table
    path trim_qc_table

  output:
    path "qc_table.txt"
    
  shell:
  '''
  awk 'BEGIN {
        OFS="\\t"
        } 
    #print header
      NR==1 {
        print $0,$0,$0,"Dropped"
        } 
    #save total read of raw to array a
      FNR>1&&NR==FNR {
        raw[FNR]=$0;
        a[FNR]=$2;
        next;
        }
    #add ratio of survived pair read
      FNR>1 && FILENAME=="qc_pair.txt" {
        pairSeq[FNR]=$2; 
        $2=$2 " (" $2/a[FNR] ")";
        b[FNR]=$0;
        next;
        } 
    #add ratio of survived unpair read
      FNR>1 && FILENAME=="qc_unpair.txt" {
        $(NF+1)=(a[FNR]-pairSeq[FNR]-$2) " (" (a[FNR]-pairSeq[FNR]-$2)/a[FNR] ")"; #add number and ratio of dropped read
        $2=$2 " (" $2/a[FNR] ")"; 
        print raw[FNR],b[FNR],$0; #join tables
        }' FS='\\t' OFS='\\t' !{raw_qc_table} !{trim_qc_table} > qc_table.txt
  
  '''
}

/*
 * Step 2a: Align reads to the reference ref_genome
 */
process READ_MAPPING_BWA {
  tag "$id"
  publishDir "$params.results/bam", mode: 'copy', \
    pattern: '*.bam*', \
    saveAs: { filename -> "${id}_$filename" }

  input: 
    path ref_genome
    val ref_genome_name
    path index
    tuple val(id), path(reads)
    val mapping_option

  output: 
    tuple \
      val(id), \
      path("markdup.bam"), \
      path("markdup.bam.bai")
    path "${id}_coverage.txt"

  """
  # read mapping
  bwa mem -R "@RG\\tID:${id}\\tSM:${id}\\tPL:Illumina" \
    ${mapping_option} $ref_genome $reads > sam

  # sort and BAM file
  samtools fixmate -O bam sam bam_fixmate
  samtools sort -O bam -o bam_sort bam_fixmate
  samtools index bam_sort
 
  # mask duplicate
  gatk --java-options "-Xmx1g" MarkDuplicates -I bam_sort \
    -O markdup.bam -M bam_markdup_metrics.txt
  
  # index BAM file
  samtools index markdup.bam

  # depth
  samtools coverage markdup.bam > coverage.txt

  sed "s/${ref_genome_name}/$id/;s/#rname/id/" coverage.txt > ${id}_coverage.txt
  
  # free up some disk space
  rm sam bam_fixmate bam_sort coverage.txt
  """
}


/*
 * Step 2b: Combine per-sample coverage info into a single file
 */
process COVERAGE_OUTPUT {
  publishDir "$params.results/bam", mode: 'copy'

  input:
    path coverage

  output:
    path "coverage.tsv"

  """
  awk '(NR == 1) || (FNR > 1)' $coverage > coverage.tsv
  """
}


/*
 * Step 3a: Call variants for each sample
 */
process CALL_VARIANTS {
  tag "$id"
  publishDir "$params.results/vcf", mode: 'copy'

  input:
    path ref_genome
    path index
    path dict
    tuple val(id), path(bam), path(bai)
    val haplotypecaller_option
 
  output: 
    path "${id}.g.vcf.gz", emit: vcf
    path "${id}.g.vcf.gz.tbi", emit: vcf_tbi

  """ 
  gatk HaplotypeCaller \
    -R $ref_genome -I $bam \
    -O ${id}.g.vcf.gz \
    -ERC GVCF \
    ${haplotypecaller_option}
  """
}


/*
 * Step 3b: Create a list of samples to allow user to exclude samples
 */
process CREATE_SAMPLE_MAP {
  publishDir "$baseDir/data", mode: 'copy'

  input:
    path vcf_all

  output:
    path "sample_map.txt", emit: sample_map

  script:
  def vcf_all_list = vcf_all.collect{ "$it\t$it" }.join('\n')
  """
  echo "${vcf_all_list}" > out1.txt
  sed 's/.g.vcf.gz\t/\t/g' out1.txt > sample_map.txt
  """
}


/*
 * Step 4: Joinly call variants from multiple samples
 */
process JOINT_GENOTYPING {
  input:
    path ref_genome
    path index
    path dict
    val ref_genome_name
    path vcf_all
    path vcf_tbi_all
    path sample_map_usr
    val genomicsdbimport_option
    val genotypegvcfs_option

  output:
    tuple path("joint.vcf.gz"), path("joint.vcf.gz.tbi")
  
  when:
    !params.stop

  script:
  def vcf_all_list = vcf_all.collect{ "-V $it" }.join(' ')
  def is_sample_map_usr = sample_map_usr.exists()

  if (is_sample_map_usr) {
    """
    gatk GenomicsDBImport \
      -R $ref_genome \
      --sample-name-map ${sample_map_usr} \
      --genomicsdb-workspace-path dbcohort \
      -L $ref_genome_name \
      ${genomicsdbimport_option}

    gatk GenotypeGVCFs -R $ref_genome -V gendb://dbcohort \
      -O joint.vcf.gz ${genotypegvcfs_option}
    """
  } else {
    """
    gatk GenomicsDBImport \
      -R $ref_genome \
      ${vcf_all_list} \
      --genomicsdb-workspace-path dbcohort \
      -L $ref_genome_name \
      ${genomicsdbimport_option}

    gatk GenotypeGVCFs -R $ref_genome -V gendb://dbcohort \
      -O joint.vcf.gz ${genotypegvcfs_option}
    """
  }
}


/*
 * Step 5: Filter variants
 */
process FILTER_VARIANTS {
  publishDir "$params.results/vcf_joint", mode: 'copy'

  input:
    path ref_genome
    path index
    path dict
    val variant_filter
    val variant_filter_name
    tuple path("joint.vcf.gz"), path("joint.vcf.gz.tbi")
    val selectvariant_option
  
  output:
    tuple path("joint_filtered.vcf.gz"), path("joint_filtered.vcf.gz.tbi")
    path "joint_filtered_stats.txt"

  when:
    !params.stop

  """
  gatk VariantFiltration \
    -R $ref_genome \
    -V joint.vcf.gz \
    -O joint_filt.vcf.gz \
    -filter "${variant_filter}" \
    --filter-name ${variant_filter_name}
  
  gatk SelectVariants \
    -V joint_filt.vcf.gz \
    -O joint_filtered.vcf.gz \
    ${selectvariant_option}

  bcftools stats -F $ref_genome joint_filtered.vcf.gz > joint_filtered_stats.txt
  """
}


/*
 * Step 6: Convert filtered vcf to multiple sequence alignment
 */
process VCF_TO_FASTA {
  publishDir "$params.results/aln", mode: 'copy'

  input:
    tuple path(vcf), path(tbi)

  output:
    path "pos.txt"
    path "aln.fasta"

  when:
    !params.stop
  
  script:
  $/
  gatk VariantsToTable -V $vcf -O gt -F POS -GF GT

  # transpose genotype table from gatk
  datamash transpose --output-delimiter=, < gt > tmp

  # create alignment file
  sed -i.tmp '1d' tmp

  sed $'s/^/>/;s/.GT,/\\\n/g;s/,//g;s/[\.\*]/-/g' tmp > aln.fasta

  # cleanup
  rm gt tmp tmp.tmp

  # optional: positions with ref variants
  gatk VariantsToTable -V $vcf -O pos.txt -F POS -F REF
  /$
}
