unzip -p !{input} !{input.baseName}/fastqc_data.txt > !{input.baseName}

  awk -F '\\t' '/Filename/ {
        flag=1
    }; 
    flag {print $2}; 
    /%GC/ {
        flag=0
    }' !{input.baseName} > !{input.baseName}.txt \

  sed -n '/Per base sequence quality/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' >> !{input.baseName}.txt \

  sed -n '/Per sequence quality scores/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR==3 {
          print $1
      }; 
      {
          max=lastline;lastline=$1
      }; 
      END {
          print max
      }' >> !{input.baseName}.txt \

  sed -n '/Per base sequence content/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2&&NF>1{
      for (i=2;i<=NF;i++){
          a[i]+=$i;
          }
      } 
      END {
          print a[2]/(NR-3),
          "\\n",a[3]/(NR-3),
          "\\n",a[4]/(NR-3),
          "\\n",a[5]/(NR-3) 
          }' >> !{input.baseName}.txt \

  sed -n '/Per base N content/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2 && NF>1{
          n++;sum+=$2
      } 
      END {
          print sum/(NR-3)
      }' >> !{input.baseName}.txt \

  sed -n '/Total Deduplicated Percentage/,/END/p' !{input.baseName} |
  awk -F '\\t' 'NR==1 {print $2}' >> !{input.baseName}.txt \

  sed -n '/Adapter Content/,/END/p' !{input.baseName} | 
  awk -F '\\t' 'NR>2&&NF>1 {
      for (i=2;i<=NF;i++){
          a[i]+=$i;
          }
      } 
      END {
          print a[2]/(NR-3),
          "\\n",a[3]/(NR-3),
          "\\n",a[4]/(NR-3),
          "\\n",a[5]/(NR-3),
          "\\n",a[6]/(NR-3)
      }' >> !{input.baseName}.txt \

  paste -s !{input.baseName}.txt > !{input.baseName}_Trimmed.txt
  
  '''