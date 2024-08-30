#!/usr/bin/bash
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-12-19, 13:59:37
# modified date: 2024-02-07, 00:08:56
shopt -s expand_aliases

if [ $# -lt 2 ];then
    echo -e "Usage: $0 <fq|fq1,fq2> <out_prefix>"
    echo -e "Example:\n ${0##*/} xx.fq.gz out_prefix"
    echo -e " ${0##*/} xx_1.fq.gz,xx_2.fq.gz out_prefix"
    exit 2
fi

fq=$1
out=$2
out_p=$(basename $2)
out_d=$(dirname $2)
klist='21,33,55,77,99'
trds=56
alias spades.py=/share/data1/software/SPAdes-3.15.5/bin/spades.py

( [ ! -d $out_d ] ) && mkdir -p $out_d

if [ -f $out.log ] && grep -q 'SPAdes pipeline finished' $out.log ;then
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] INFO Existing and skip sample: $out_p."
    exit 0
elif [ -d $out ] && ! grep -q 'SPAdes pipeline finished' $out/spades.log ;then
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] TASK Removing the folder: $out_p."
    rm -rf $out
fi

if [[ $fq =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    spades.py --isolate -1 $fq1 -2 $fq2 -k $klist -o $out -t $trds > /dev/null 2>/dev/null || { rm -r $out && exit 1; }
else
    spades.py --isolate -s $fq -k $klist -o $out -t $trds > /dev/null 2>/dev/null || { rm -r $out && exit 1; }
fi
mv $out/scaffolds.fasta $out.fa
mv $out/spades.log $out.log
chmod 444 $out.fa $out.log
rm -rf $out
