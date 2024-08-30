#!/usr/bin/bash
shopt -s expand_aliases

( [ $# -lt 3 ] ) && echo -e "Usage: bash $0 [contigs_fasta] [min_length] [out_prefix] [threads:80]" && exit -1

fa=$1; len=$2; out_p=$3; trds=${4:-80}
alias seqkit=/share/data1/software/binary/seqkit
alias hmmsearch=/share/data1/software/miniconda3/envs/checkv/bin/hmmsearch
alias diff_head=/share/data1/mjx/bin/diff_head
alias flow_prodigal.sh=/share/data1/mjx/bin/flow_prodigal.sh
alias filter_hmm2.pl=/share/data1/mjx/bin/filter_hmm2.pl
alias parallel=/share/data1/software/parallel-20230822/bin/parallel
ENV_PATH_CHECKV=/share/data1/software/miniconda3/envs/checkv
ENV_PATH_VIBRANT=/share/data1/software/miniconda3/envs/VIBRANT-1.2.1
ENV_PATH_DEEPVIRFINDER=/share/data1/software/miniconda3/envs/deepvirfinder
DVF_MODEL=/share/data1/software/miniconda3/DeepVirFinder/models
BUSCO_HMMS=/share/data1/database/busco_v5/busco_downloads/lineages/bacteria_odb10/bacteria_odb10.hmms
BUSCO_CUTOFF=/share/data1/database/busco_v5/busco_downloads/lineages/bacteria_odb10/scores_cutoff
export CHECKVDB=/share/data1/database/checkv-db-v1.5/
export VIBRANT_DATA_PATH=/share/data1/database/vibrant/

print_log(){ 
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] $*" 
}

( [ ! -d $out_p.phages ] ) && mkdir -p $out_p.phages
cd $out_p.phages
fa=../$fa
( [ ! -d tmp ] ) && mkdir -p tmp

# filter out sequences that are unlikely to be bacteriophage firstly.
( [ ! -d 1_checkv ] ) && mkdir 1_checkv && \
    print_log TASK - Step [1/6]: filter out sequences that are unlikely to be bacteriophage by using checkv firstly. > $out_p.log && \
    source activate $ENV_PATH_CHECKV >/dev/null 2>/dev/null && \
    checkv end_to_end $fa 1_checkv -t $trds >> $out_p.log 2>&1 && \
    awk '($3=="No"||$1=="contig_id") && !($7/($5+1e-16)>0.5 && ($6==0 || ($6>0 && $7/($6+1e-16)>5)))' 1_checkv/quality_summary.tsv > tmp/quality_summary.tsv && \
    n=`awk 'END{printf "%d\n",NR/10+1}' tmp/quality_summary.tsv` && \
    cut -f1 tmp/quality_summary.tsv | split -d -l $n - tmp/list.part_ && \
    find tmp/list.part_* | parallel -j 10 seqkit grep --quiet -f {} -o {}.fa 1_checkv/viruses.fna && \
    cat tmp/list.part_*.fa > tmp/virus.fa && rm -r tmp/list.part_* && \
    sed 's/ .*//' 1_checkv/proviruses.fna >> tmp/virus.fa && \
    mkdir -p 1_checkv/recheckv && checkv end_to_end 1_checkv/proviruses.fna 1_checkv/recheckv -t $trds >> $out_p.log 2>&1 && \
    awk 'NR!=1 && $2>'$len'' 1_checkv/recheckv/quality_summary.tsv >> tmp/quality_summary.tsv && \
    cat 1_checkv/tmp/proteins.faa 1_checkv/recheckv/tmp/proteins.faa > tmp/virus.faa && \
    conda deactivate || { rm -rf 1_checkv && exit 2; }
    
    # 初步过滤基本不可能是病毒的序列，减少之后病毒序列预测的数据量
    # 不要宿主基因数量占总基因数量一半以上的序列
    # 不要病毒基因为0的序列，或者是有病毒基因，但是宿主基因数量是病毒基因5倍的序列
    # 前病毒的序列补充到过滤的序列中，这些原病毒序列已经被剪切掉宿主序列了
    # 针对前病毒重新跑一下checkv，把病毒序列从前病毒序列中提取出来，也把质量报告合并到之前的指控报告中

# vibrant predict bacteriophage.
( [ ! -d 2_VIBRANT ] ) && mkdir 2_VIBRANT && \
    print_log TASK - Step [2/6]: vibrant predict bacteriophage. >> $out_p.log && \
    source activate $ENV_PATH_VIBRANT >/dev/null 2>/dev/null && \
    VIBRANT_run.py -l $len -f nucl -no_plot -t $trds -i tmp/virus.fa -folder 2_VIBRANT >> $out_p.log 2>&1 && \
    conda deactivate || { rm -rf 2_VIBRANT && exit 2; }

# deepvirfinder predict bacteriophage.
( [ ! -d 3_dvf ] ) && mkdir 3_dvf && \
    print_log TASK - Step [3/6]: deepvirfinder predict bacteriophage. >> $out_p.log && \
    source activate $ENV_PATH_DEEPVIRFINDER >/dev/null 2>/dev/null && \
    dvf.py -i tmp/virus.fa -o 3_dvf -l $len -m $DVF_MODEL -c 8 >> $out_p.log 2>&1 && \
    conda deactivate || { rm -rf 3_dvf && exit 2; }

## Merge result.
( [ ! -f tmp/result.faa ] ) && print_log TASK - Step [4/6]: combine result. >> $out_p.log && \
    awk '$6>$7 && $8!="Not-determined"' tmp/quality_summary.tsv > tmp/checkv_quality_summary.tsv && \
    perl -ne 's/_fragment.*//;print $_;' 2_VIBRANT/VIBRANT_virus/VIBRANT_phages_virus/virus.phages_combined.txt |\
    sort --parallel=5 -u > tmp/VIBRANT_virus.tsv && \
    awk '{if($3>0.9 && $4<0.01) print $1}' 3_dvf/virus.fa_gt${len}bp_dvfpred.txt > tmp/dvf_pred.tsv && \
    cat tmp/checkv_quality_summary.tsv tmp/VIBRANT_virus.tsv tmp/dvf_pred.tsv > tmp/combine_predicted_result.tsv && \
    diff_head tmp/combine_predicted_result.tsv 1 tmp/checkv_quality_summary.tsv 1 tmp/result.tsv 2>/dev/null && \
    cut -f1 tmp/result.tsv | sed '1d;s/$/_/' | seqkit grep --quiet -r -f - tmp/virus.faa > tmp/result.faa || \
    { ! rm -rf tmp/checkv_quality_summary.tsv tmp/VIBRANT_virus.tsv 5_gather/dvf_pred.tsv tmp/combine_predicted_result.tsv tmp/result.tsv tmp/result.faa && exit 2; }

## decontaminated by BUSCO
( [ ! -d 4_BUSCO ] ) && mkdir 4_BUSCO && \
    print_log TASK - Step [5/6]: decontaminated by BUSCO. >> $out_p.log && \
    hmmsearch --tblout 4_BUSCO/hmm.out --cpu $trds --noali $BUSCO_HMMS tmp/result.faa >/dev/null 2>/dev/null && \
    filter_hmm2.pl 4_BUSCO/hmm.out $BUSCO_CUTOFF 4_BUSCO/hmm.out.cut tmp/result.faa 4_BUSCO/hmm.out.stat >/dev/null 2>&1 && \
    awk '$3/$2<0.05' 4_BUSCO/hmm.out.stat > 4_BUSCO/result.tsv || \
    { ! rm -rf 4_BUSCO && exit 2; }

## Result
( [ ! -f $out_p.virus.fa ] ) && \
    print_log TASK - Step [6/6]: output bacteriophage quality summary and sequences. >> $out_p.log && \
    diff_head 4_BUSCO/result.tsv 1 tmp/checkv_quality_summary.tsv 1 $out_p.virus.qs 2>/dev/null && \
    sed '1d' $out_p.virus.qs | cut -f1 | seqkit grep --quiet -f - tmp/virus.fa -o $out_p.virus.fa || \
    { rm $out_p.virus.qs $out_p.virus.fa && exit 2; }
