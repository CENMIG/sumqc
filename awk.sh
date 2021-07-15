extract_fastqc () {
    if [ -z $1  ];
    then
        echo "Supplies one or more fastqc result in zip file."
    fi
    unset FileName BaseNameSignal OmitSignal OPTIND
    FileName="FileName\t"
    while getopts "bo" opt;
    do
        case $opt in
            #Print filename as in basename
            b) BaseNameSignal='1';;

            #Omit filename column
            o) OmitSignal='1'
               FileName=""
               ;;
        esac
    done
    shift "$(( OPTIND - 1 ))"

    echo -e "${FileName}TotalSeq\tLength\tGC\tAvgQScore (min,max)";
    for file in $@
    do
        trimfile=${file##*/}
        unzip -p $file ${trimfile%.*}/fastqc_data.txt | \
        awk -F '\t'  -v omit=$OmitSignal -v base=$BaseNameSignal -v name=${trimfile%.*} \
        '/Filename/{if ( length(base)==0 ) {file=$2} else {gsub(/_[12].*/,"");file=$2}}
        /Total Sequences/{seq=$2}
        /Sequence length/{len=$2}
        /%GC/{gc=$2}
        /Per sequence quality scores/,/END/{nr[NR]=$1; n++ ;nrr[n]=$1 ; min=nrr[3]; max=nr[NR-1] ;$3=$1*$2; sumn+=$2;sum+=$3 ;}
    END{if (length(omit)==0){
        if ( sumn==0 ) {printf "%s\t%s\t%s\t%.3f\t%.3f (%s,%s)\n",file,0,0,0,0,0,0 ; exit};
        printf "%s\t%s\t%s\t%.3f\t%.3f (%s,%s)\n",file,seq,len,gc,sum/sumn,min,max
         }
    else {
        if ( sumn==0 ) {printf "%s\t%s\t%.3f\t%.3f (%s,%s)\n",0,0,0,0,0,0 ; exit};
        printf "%s\t%s\t%.3f\t%.3f (%s,%s)\n",seq,len,gc,sum/sumn,min,max
         }}' OFS='\t'
    done
}
#print header then continue with whatever command
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
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
