#!/usr/bin/bash
shopt -s expand_aliases

if [ $# -lt 3 ];then
    echo "$0 [fq|fq1,fq2] [kk2_db] [out_prefix] [read_len:-150]"
    echo -e " \033[31mOptional k2_db:\033[0m:"
    echo -e "  k2_standard_16gb [refeq archaea, bacteria, viral, plasmid, human, UniVec_Core]"
    echo -e "  k2_pluspf [standard plus Refeq protozoa & fungi]"
    echo -e "  k2_16S_SILVA138 [16s sliva RNA indexes]"
    echo -e "  k2_uhgg_v2.0.2 [UHGG]"
    exit 2
fi

fq=$1
db=$2
p=$3
len=${4:-150}
trds=16
alias kraken2=/share/data1/software/miniconda3/envs/kraken2/bin/kraken2
alias bracken=/share/data1/software/miniconda3/envs/kraken2/bin/bracken

# index
declare -A k2_db
k2_db["k2_standard_16gb"]="/share/data1/database/kraken2_db/k2_standard_16gb_20230605"
k2_db["k2_pluspf"]="/share/data1/database/kraken2_db/k2_pluspf_20230605"
k2_db["k2_16S_SILVA138"]="/share/data1/database/kraken2_db/16S_SILVA138_k2db_20200326"
k2_db["k2_uhgg_v2.0.2"]="/share/data1/database/kraken2_db/k2_uhgg_v2.0.2_20240129"
db=${k2_db["${db}"]}

if [[ $fq =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    kraken2 --db $db --output $p.out --report $p.kk2 --threads $trds --use-names --gzip-compressed $fq1 $fq2 >> $p.log 2>&1
else
    kraken2 --db $db --output $p.out --report $p.kk2 --threads $trds --use-names --gzip-compressed $fq >> $p.log 2>&1
fi

bracken -l P -r $len -d $db -i $p.kk2 -o $p.bk.p -w $p.bk.p.out >> $p.log 2>&1 && rm $p.bk.p.out 
# bracken -l F -r $len -d $db -i $p.kk2 -o $p.bk.f -w $p.bk.f.out >> $p.log 2>&1 && rm $p.bk.f.out
bracken -l G -r $len -d $db -i $p.kk2 -o $p.bk.g -w $p.bk.g.out >> $p.log 2>&1 && rm $p.bk.g.out
bracken -l S -r $len -d $db -i $p.kk2 -o $p.bk.s -w $p.bk.s.out >> $p.log 2>&1 && rm $p.bk.s.out
rm $p.out 

