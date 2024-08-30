#!/usr/bin/bash
shopt -s expand_aliases

if [ $# -ne 3 ];then
    echo -e "$0 <fq|fq1,fq2> <host> <out_prefix>\e[0m"
    echo -e " \033[31mOptional host:\033[0m:"
    printf "%16s:    %-20s %-20s\n" human Homo_sapiens GCF_000001405.40_GRCh38
    printf "%16s:    %-20s %-20s\n" pig Sus_scrofa GCA_000003025.6
    printf "%16s:    %-20s %-20s\n" chicken Gallus_gallus GCF_016699485.2
    exit 2
fi

fq=$1
host=$2
p=$3
trds=16
alias fastp="/share/data1/software/binary/fastp"
alias bowtie2="/share/data1/software/bowtie2-2.5.1/bowtie2"
# index
declare -A host_index
host_index["human"]="/share/data1/database/genome_host/Homo_sapiens_human/Homo_sapiens_human"
host_index["pig"]="/share/data1/database/genome_host/Sus_scrofa_pig/Sus_scrofa_pig"
host_index["chicken"]="/share/data1/database/genome_host/Gallus_gallus_chicken/Gallus_gallus_chicken"

idx=${host_index["${host}"]}

if [ -f ${p}_map2host.log ] && grep -q 'overall' ${p}_map2host.log;then
    echo -e "[$(date +%Y-%m-%d\ %H:%M:%S)] INFO Existing and skip sample: $p."
    exit 0
fi
# fastp and remove host
if [[ ${fq} =~ "," ]];then
    fq1=$(echo $fq | cut -d "," -f1)
    fq2=$(echo $fq | cut -d "," -f2)
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g \
        -i $fq1 -I $fq2 -o ${p}_clean_1.fq.gz -O ${p}_clean_2.fq.gz -h /dev/null -j /dev/null 1>${p}_fastp.log 2>&1
    bowtie2 --end-to-end --mm --fast -p $trds -x $idx --no-head -1 ${p}_clean_1.fq.gz -2 ${p}_clean_2.fq.gz 2> ${p}_map2host.log |\
        perl -ne 'chomp;@s=split /\s+/;if($s[1]==77){print "\@$s[0]/1\n$s[9]\n+\n$s[10]\n";}elsif($s[1]==141){print STDERR "\@$s[0]/2\n$s[9]\n+\n$s[10]\n";}' \
        > >(pigz> ${p}_rmhost_1.fq.gz) 2> >(pigz > ${p}_rmhost_2.fq.gz)
    chmod 444 ${p}_fastp.log ${p}_map2host.log ${p}_rmhost_1.fq.gz ${p}_rmhost_2.fq.gz
    rm ${p}_clean_1.fq.gz ${p}_clean_2.fq.gz
else
    fastp -w $trds -q 20 -u 30 -n 5 -y -Y 30 -l 80 --trim_poly_g \
        -i $fq -o ${p}_clean.fq.gz -h /dev/null -j /dev/null 1>${p}_fastp.log 2>&1
    bowtie2 --end-to-end --mm --fast -p $trds -x $idx --no-head -U ${p}_clean.fq.gz 2> ${p}_map2host.log |\
        perl -ne 'chomp;@s=split /\s+/;if($s[1]==4){print "\@$s[0]\n$s[9]\n+\n$s[10]\n";}' \
        > >(pigz > ${p}_rmhost.fq.gz)
    chmod 444 ${p}_fastp.log ${p}_map2host.log ${p}_rmhost.fq.gz
    rm ${p}_clean.fq.gz
fi

