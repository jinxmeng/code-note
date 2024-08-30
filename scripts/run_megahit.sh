#!/usr/bin/bash
# created date: 2023-01-01, 11:34:27
# modified date: 2024-02-11, 19:44:58

shopt -s expand_aliases

if [ $# -lt 2 ];then
    echo -e "$0 [fq|fq1,fq2] [out_prefix]"
    echo -e "Example:\n ${0##*/} xx.fq.gz out_prefix"
    echo -e " ${0##*/} xx_1.fq.gz,xx_2.fq.gz out_prefix"
    exit 2
fi

fq=$1
out=$2
out_p=$(basename $2)
out_d=$(dirname $2)
klist='21,41,61,81,101,121,141'
trds=32
alias megahit=/share/data1/software/binary/megahit

( [ ! -d $out_d ] ) && mkdir -p $out_d

if [ -f $out.log ] && grep -q 'ALL DONE' $out.log;then
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] INFO Existing and skipping sample: $out_p."
    exit 0
elif [ -d $out ] && ! grep -q 'ALL DONE' $out/log;then
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] TASK Remove the folder: $out_p."
    rm -rf $out
fi

if [[ $fq =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    megahit --k-list $klist -t $trds -1 $fq1 -2 $fq2 --out-dir $out > /dev/null 2>/dev/null || { rm -rf $out && exit 1; }
else
    megahit --k-list $klist -t $trds -r $fq --out-dir $out > /dev/null 2>/dev/null || { rm -rf $out && exit 1; }
fi
mv $out/final.contigs.fa $out.fa
mv $out/log $out.log
chmod 444 $out.fa $out.log
rm -rf $out
