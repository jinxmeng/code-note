#!/usr/bin/bash
shopt -s expand_aliases

( [ $# -ne 2 ] ) && { echo -e "$0 [SRR...] [out_directory]\nFor example: $0 SRR16020678 out_directory" && exit 2; }

srr=$1; out=$2
( [[ ! $srr == *"SRR"* ]] ) && { echo -e "Skip $srr .." && exit 0; }
( `ls $out | grep -q ${srr}*gz` ) && { echo -e "Skip sample: $srr .."; exit 0; }

alias prefetch=/share/data1/software/sratoolkit.3.0.2/bin/prefetch
alias vdb-validate=/share/data1/software/sratoolkit.3.0.2/bin/vdb-validate
alias fastq-dump=/share/data1/software/sratoolkit.3.0.2/bin/fastq-dump

[ ! -d $out ] && mkdir $out
prefetch $srr -X 200G -C yes --output-file $out/$srr.sra
vdb-validate $out/$srr.sra
fastq-dump $out/$srr.sra --split-3 --gzip -O $out
rm $out/$srr.sra
