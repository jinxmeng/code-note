#!/usr/bin/bash
# encoding: utf-8
# author: Jinxin Meng
# created date: 2024-02-24, 16:21:04
# modified date: 2024-05-14, 21:32:18
shopt -s expand_aliases 

# 2024-05-14: add -a parameter.

usage() {
    echo -e "Usage: bash $0 [-i fq|fq1,fq2 or *.gz format] [-r contigs.fa] [-t minimap2] [-p 16] [-o out_file]"
    echo -e "  Option:"
    echo -e "  -t CHR    specify a tool for short-reads alignment, including bwa, bwa-mem2, bowtie2 and minimap2 [minimap2]"
    echo -e "  -p NUM    threads number [16]"
    echo -e "  -i CHR    short-reads files, fq for single-end or fq1,fq2 for paired-end"
    echo -e "  -r CHR    reference fasta, such as contigs.fa -r/-x conflicting"
    echo -e "  -x CHR    reference index, -r/-x conflicting, format following:"
    echo -e "            *.mmi for minimap2, need whole file name"
    echo -e "            *.amb, *.ann, *.bwt, *.pac and *.sa for bwa, need suffix"
    echo -e "            *.0123, *.amb, *.ann, *.bwt.2bit.64 and *.pac for bwa-mem2, need suffix"
    echo -e "            *.1.bt2, *.2.bt2, *.3.bt2, *.4.bt2, *.rev.1.bt2, *.rev.2.bt2 for bowtie2, need suffix"
    echo -e "  -a CHR    append parameter to alignment tools"
    echo -e "  -o CHR    output SAM file prefix"
    echo -e "  -b        output sort BAM file"
    echo -e "  -c        clean tmp file, such as index file for bwa, bwa-mem2, bowtie2"
    echo -e "  -h        print manual"
    echo -e "  Example:"
    echo -e "  $0 -i xx.fq.gz -r xx.fa -o xx -t minimap2 -p 16"
    echo -e "  $0 -i xx_1.fq.gz,xx_2.fq.gz -x xx -o xx -t bwa-mem2 -p 24 -b -c"
    echo -e "  Other:"
    echo -e "  samtools used to process sam or bam file [samtools, sambamba]"
    exit -1
}

print_log(){
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] $*"
}

[ $# -eq 0 ] && { usage && exit -1; }

trds=16; tool=minimap2; bam="F"; clean="F"; add="";
while getopts "t:p:i:r:x:m:a:o:bch" arg; do
    case $arg in
        t) tool="$OPTARG";;
        p) trds="$OPTARG";;
        i) fq="$OPTARG";;
        r) fa="$OPTARG";;
        x) index="$OPTARG";;
        a) add="$OPTARG";;
        o) out="$OPTARG";;
        b) bam="T";;
        c) clean="T";;
        h) { usage && exit -1; };;
        ?) { usage && exit -1; };;
    esac
done

( [ ! $fq ] || ([ ! $fa ] && [ ! $index ]) || [ ! $out ] ) && { print_log fq, out and fa/index need to be given .. && exit -1; } # parameters must have fq, fa, out
[ -f $out.log ] && grep -q 'Program success' $out.log && { print_log INFO skipping sample: $out .. && exit 0; } # skip run
([ $index ] && [ $fa ]) && unset $fa # if both index and fa be given, then unset fa
[[ $fq =~ "," ]] && fq1=$(echo $fq | cut -d "," -f1); fq2=$(echo $fq | cut -d "," -f2) # split fq

alias bwa=/share/data1/software/bwa-0.7.17/bwa
alias bwa-mem2=/share/data1/software/bwa-mem2-2.2.1/bwa-mem2
alias minimap2=/share/data1/software/minimap2-2.26/minimap2
alias bowtie2-build=/share/data1/software/bowtie2-2.5.1/bowtie2-build
alias bowtie2=/share/data1/software/bowtie2-2.5.1/bowtie2
alias samtools=/share/data1/software/samtools/bin/samtools
alias sambamba=/share/data1/software/binary/sambamba

# commend
cmd_log="$out.log"
cmd_bam="| samtools view -@ $trds -F 3584 -bS | samtools sort -@ $trds -o $out.sort.bam 2>/dev/null"
# cmd_bam="| sambamba view -q -S -t $trds -F \"not (failed_quality_control or duplicate or supplementary)\" -f bam -o /dev/stdout /dev/stdin | sambamba sort -q -t $trds -o $out.sort.bam /dev/stdin"
cmd_index_existing="print_log TASK $tool index existing .. >>$cmd_log"
cmd_indexing="print_log TASK $tool indexing .. >>$cmd_log"
cmd_mapping="print_log TASK $tool mapping .. >>$cmd_log"

if [ $tool == "minimap2" ];then
    [ $index ] && cmd="$cmd_index_existing && $cmd_mapping"
    [ ! $index ] && { cmd="$cmd_indexing && minimap2 -d ${fa%.*}.mmi $fa -I30g 2>>$cmd_log && $cmd_mapping" && index=${fa%.*}.mmi; }
    cmd="$cmd && minimap2 -t $trds -ax sr $index $add"
    index_list=$index
    [[ $fq =~ "," ]] && cmd="$cmd $fq1 $fq2 2>>$cmd_log" || cmd="$cmd $fq 2>>$cmd_log"
    [ $bam == "T" ] && cmd="$cmd $cmd_bam" || cmd="$cmd -o $out.sam"
    [ $clean == "T" ] && cmd="$cmd && rm $index_list 2>/dev/null"
fi

if [ $tool == "bwa" ];then
    [ $index ] && cmd="$cmd_index_existing && $cmd_mapping"
    [ ! $index ] && { cmd="$cmd_indexing && bwa index $fa -p ${fa%.*} 2>>$cmd_log && $cmd_mapping" && index=${fa%.*}; }
    cmd="$cmd && bwa mem -t $trds $index $add"
    index_list="$index.amb $index.ann $index.bwt $index.pac $index.sa"
    [[ $fq =~ "," ]] && cmd="$cmd $fq1 $fq2 2>>$cmd_log" || cmd="$cmd $fq 2>>$cmd_log"
    [ $bam == "T" ] && cmd="$cmd $cmd_bam" || cmd="$cmd -o $out.sam"
    [ $clean == "T" ] && cmd="$cmd && rm $index_list 2>/dev/null"
fi

if [ $tool == "bwa-mem2" ];then
    [ $index ] && cmd="$cmd_index_existing && $cmd_mapping"
    [ ! $index ] && { cmd="$cmd_indexing && bwa-mem2 index $fa -p ${fa%.*} >>$cmd_log 2>>$cmd_log && $cmd_mapping" && index=${fa%.*}; }
    cmd="$cmd && bwa-mem2 mem -t $trds $index $add"
    index_list="$index.0123 $index.amb $index.ann $index.bwt.2bit.64 $index.pac"
    [[ $fq =~ "," ]] && cmd="$cmd $fq1 $fq2 2>>$cmd_log" || cmd="$cmd $fq 2>>$cmd_log"
    [ $bam == "T" ] && cmd="$cmd $cmd_bam" || cmd="$cmd -o $out.sam"
    [ $clean == "T" ] && cmd="$cmd && rm $index_list 2>/dev/null"
fi

if [ $tool == "bowtie2" ];then
    [ $index ] && cmd="$cmd_index_existing && $cmd_mapping"
    [ ! $index ] && { cmd="$cmd_indexing && bowtie2-build $fa ${fa%.*} --threads $trds --quiet >>$cmd_log 2>>$cmd_log && $cmd_mapping" && index=${fa%.*}; }
    cmd="$cmd && bowtie2 --end-to-end --mm -p $trds -x $index $add"
    index_list="$index.1.bt2 $index.2.bt2 $index.3.bt2 $index.4.bt2 $index.rev.1.bt2 $index.rev.2.bt2"
    [[ $fq =~ "," ]] && cmd="$cmd -1 $fq1 -2 $fq2 2>>$cmd_log" || cmd="$cmd -U $fq 2>>$cmd_log"
    [ $bam == "T" ] && cmd="$cmd $cmd_bam" || cmd="$cmd -S $out.sam"
    [ $clean == "T" ] && cmd="$cmd && rm $index_list 2>/dev/null"
fi
print_log CMD $cmd >$out.log
eval $cmd || { rm -rf $out.log $index_list 2>/dev/null && exit -1; }
print_log INFO Program success .. >>$out.log
