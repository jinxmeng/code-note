#!/usr/bin/bash

if [ $# -lt 2 ] || [[ $1 =~ "-h" ]];then
    echo "$0 <fq>|<fq1,fq2> <out_prefix> [nreads:999999999]"
    echo "  MetaPhlAn version 4.0.6 (1 Mar 2023)"
    echo "  Database in: /share/data1/database/mpa4/mpa_vOct22_CHOCOPhlAnSGB_202212"
    exit 127
fi

set -f
set -e # 如果出错，就不再执行下一步

fq1=$1 # 哪一端的reads无所谓
nreads=${3:-999999999}

if [[ $1 =~ "," ]];then
    fq1=`echo $1|cut -d "," -f 1`
    fq2=`echo $1|cut -d "," -f 2`
    if [ ! -f $fq1 ] || [ ! -f $fq2 ];then
        echo "No such file or directory"
        exit 127
    fi
fi

if [ $# -eq 1 ]||[ ! -f $fq1 ];then
    echo "$0 <fq|fq1,fq2> <output_prefix>"
    exit 127
fi

db_path=/share/data1/database/mpa4/mpa_vOct22_CHOCOPhlAnSGB_202212
sotf_bw=/share/data1/software/miniconda3/envs/metaphlan/bin/bowtie2
soft_mp=/share/data1/software/miniconda3/envs/metaphlan/bin/metaphlan

out_f=$2

${sotf_bw} --fast --mm -u $nreads -p 16 --sam-no-hd --sam-no-sq --no-unal -S ${out_f}.sam -x ${db_path} -U ${fq1} 2> ${out_f}.bowtie.log \
    && ${soft_mp} ${out_f}.sam --input_type sam --nreads $nreads -o ${out_f}.tsv \
    && rm ${out_f}.sam
