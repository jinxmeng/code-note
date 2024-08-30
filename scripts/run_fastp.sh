#!/usr/bin/bash
shopt -s expand_aliases

( [ $# -ne 2 ] ) &&  { echo -e "Usage: $0 [fq|fq1,fq2] [out_prefix]" && exit 2; }

fq=$1
out=$2
trds=16
alias fastp="/share/data1/software/binary/fastp"

( [ -e ${out}_fastp.log ] ) && { echo -e "Skip sample: ${out##*/} .." && exit 0; }

if [[ $fq =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g -i $fq1 -I $fq2 -o ${out}_clean_1.fq.gz -O ${out}_clean_2.fq.gz \
        -h /dev/null -j /dev/dell 1>${out}_fastp.log 2>${out}_fastp.log || { rm ${out}_* && exit 2; }
else
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g -i $fq -o ${out}_clean.fq.gz \
        -h /dev/null -j /dev/null 1>${out}_fastp.log 2>${out}_fastp.log || { rm ${out}_* && exit 2; }
fi
chmod 444 ${out}_*
