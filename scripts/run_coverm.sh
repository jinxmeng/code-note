#!/usr/bin/bash
shopt -s expand_aliases

[ $# -lt 3 ] && { echo "Usage: bash $0 [*sort.bam] [genome_list] [out_prefix] [covered_fraction: 0 <0-100>]" && exit 127; }

bam=$1; fa=$2; out=$3; cf=${4:-0}
alias coverm=/share/data1/software/binary/coverm

[ -f $out.cvm ] && { echo "skipping sample: $out .." && exit 0; }

if [ $cf -eq 0 ];then
    coverm genome --bam-files $bam --genome-fasta-list $fa --methods length count covered_fraction covered_bases tpm rpkm relative_abundance \
        --min-covered-fraction 0 --output-file $out.cvm || { rm $out.cvm && exit -1; }
else
    coverm genome --bam-files $bam --genome-fasta-list $fa --methods covered_fraction covered_bases tpm rpkm relative_abundance \
        --min-covered-fraction $cf --output-file $out.cvm || { rm $out.cvm && exit -1; }
fi
