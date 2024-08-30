#!/usr/bin/bash
shopt -s expand_aliases

( [ $# -lt 3 ] ) && { echo -e "Usage: bash $0 [contigs_fasta] [min_length] [out_prefix] [threads:112]" && exit 2; }

fa=$1
len=$2
out=$3
trds=${4:-112}
ENV_PATH_CHECKV=/share/data1/software/miniconda3/envs/checkv
ENV_PATH_VIRSORTER2=/share/data1/software/miniconda3/envs/vs2

( [ ! -d $out.phages ] ) && mkdir -p $out.phages
cd $out.phages && fa=../$fa
# ( [ ! -d tmp ] ) && mkdir -p tmp

if [ ! -d vs2.out ];then
    source activate $ENV_PATH_VIRSORTER2 >/dev/null 2>/dev/null && \
        virsorter run --keep-original-seq -w vs2.out -i $fa --include-groups dsDNAphage,ssDNA --min-score 0.5 --min-length $len -j $trds all > $out.log 2>&1 && \
        conda deactivate || \
        { rm -rf vs2.out >/dev/null 2>/dev/null && exit -1; }
        # seqkit fx2tab --quiet -n -l vs2.out/final-viral-combined.fa | \
        # perl -e '%l;%n;while(<>){chomp;@s=split/\t/;$_=~/(\S+)\|\|/;$x=$1;if(!exists $l{$x}){$l{$x}=$s[1];$n{$x}=$s[0]}else{if($s[1]>$l{$x}){$l{$x}=$s[1];$n{$x}=$s[0]}}} for (keys %n){print "$_\t$n{$_}\n"}' > tmp/vs2.out.seqid && \
        # cut -f2 tmp/vs2.out.seqid | seqkit grep --quiet -f - vs2.out/final-viral-combined.fa | seqkit replace --quiet -p "\|\|.*" -r "" -o tmp/vs2.out.fa || \
        # { rm -rf vs2.out tmp/vs2.out.seqid tmp/vs2.out.fa >/dev/null 2>/dev/null && exit -1; }
fi
if [ ! -d ckv.out ];then
    source activate $ENV_PATH_CHECKV >/dev/null 2>/dev/null && \
        checkv end_to_end -t $trds vs2.out/final-viral-combined.fa ckv.out >> $out.log 2>&1 && \
        conda deactivate && \
        perl -lne 'if(/(>\S+)/){print "$1"}else{print "$_"}' ckv.out/proviruses.fna | \
        cat ckv.out/viruses.fna - > virus.fa || \
        { rm -rf ckv.out virus.fa >/dev/null 2>/dev/null && exit -1; }
fi
echo "[$(date +%Y-%m-%d\ %H:%M:%S)] INFO: bacteriophages discovery success .." >> $out.log
