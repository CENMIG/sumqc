extract_fastqc () {
    echo -e "FileName\tTotalSeq\tPoorQualflag\tLength\tGC\tAvgQScore (min,max)";
    for file in $@
    do
        unzip -p $file ${file%.*}/fastqc_data.txt | \
        awk -F '\t' '/Filename/{file=$2}
        /Total Sequences/{seq=$2}
        /Sequences flagged as poor quality/{flg=$2}
        /Sequence length/{len=$2}
        /%GC/{gc=$2}
        /Per sequence quality scores/,/END/{nr[NR]=$1; n++ ;nrr[n]=$1 ; min=nrr[3]; max=nr[NR-1] ;$3=$1*$2; sumn+=$2;sum+=$3 ;}
        END{printf "%s\t%s\t%s\t%s\t%.3f\t%.3f (%s,%s)\n",file,seq,flg,len,gc,sum/sumn,min,max}';
    done
    }

calculate_percent () {
    #hard code col number $2-raw $8-pair $14-unpair $19-droppped (add)
    paste $@ | \
        awk -F '\t' 'NR==1{$19="dropped"}
        NR>1{pratio=($2-$8)*100/$2;
        uratio=($2-$14)/$2*100;
        drop=$2-($8-$14);
        dratio=($2-$8-$14)/$2*100;
        format="%d (%.2f)";
        $8=sprintf(format,$8,pratio);
        $14=sprintf(format,$14,uratio);
        $19=sprintf(format,drop,dratio);}
        {print}' OFS='\t';

}
