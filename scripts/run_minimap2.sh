#!/usr/bin/bash
#Keep identity >95%

( [ $# -ne 3 ] ) &&  { echo -e "Usage: $0 [fq|fq1,fq2] [*.mmi] [out_prefix]" && exit 2; }
( [ -e $3.sort.bam ] ) && { echo -e "Skip sample: ${3##*/} .." && exit 0; }

if [[ $1 =~ "," ]];then
    fq1=$(echo $1 | cut -d "," -f1)
    fq2=$(echo $1 | cut -d "," -f2)
    minimap2 --MD -t 24 -ax sr $2 $fq1 $fq2 2>>$3.log |\
        perl -e 'while(<>){if(/^@/){print "$_"; next};chomp;@l=split/\t/;;next if ($l[1] & 0x4) != 0;$ref_length=0;$match_counts=0;$cov_length=0;while($l[5]=~/(\d+)[M=XID]/g){$cov_length+=$1};while($l[5]=~/(\d+)[MDN=X]/g){$ref_length+=$1};foreach $k(@l[11..$#l]){if($k=~/MD:Z:/){while($k=~/(\d+)/g){$match_counts+=$1}}};$identity=$match_counts/$cov_length*100;print "$_\n" if $identity>=95}' |\
        samtools view -@ 24 -bS | samtools sort -@ 24 -o $3.sort.bam -
else
    minimap2 --MD -t 24 -ax sr $2 $1 2>>$3.log |\
        perl -e 'while(<>){if(/^@/){print "$_"; next};chomp;@l=split/\t/;;next if ($l[1] & 0x4) != 0;$ref_length=0;$match_counts=0;$cov_length=0;while($l[5]=~/(\d+)[M=XID]/g){$cov_length+=$1};while($l[5]=~/(\d+)[MDN=X]/g){$ref_length+=$1};foreach $k(@l[11..$#l]){if($k=~/MD:Z:/){while($k=~/(\d+)/g){$match_counts+=$1}}};$identity=$match_counts/$cov_length*100;print "$_\n" if $identity>=95}' |\
        samtools view -@ 24 -bS | samtools sort -@ 24 -o $3.sort.bam -
fi
