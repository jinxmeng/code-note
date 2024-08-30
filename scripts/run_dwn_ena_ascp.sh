#!/usr/bin/bash
shopt -s expand_aliases

( [ $# -ne 2 ]) && { echo -e "$0 [filereport_read_run_xxx_tsv.txt] [out_directory]\n Notice: dwn_list must contain \"fastq_aspera\"" && exit 2; }

export in_f=$1
out=$2

alias parallel=/share/data1/software/parallel-20230822/bin/parallel
alias ascp=/home/mjx/.aspera/connect/bin/ascp
openssh=/home/mjx/.aspera/connect/etc/asperaweb_id_dsa.openssh

/usr/bin/perl << 'EOF' > tmp_dwn_list
use strict;
use warnings;
open I, "<$ENV{in_f}"; # perl 获取 shell变量
my $col;
while (<I>) {
    chomp;
    my @s = split /\t/;    
    if (/fastq_aspera/) {
        ($col) = grep { $s[$_] eq 'fastq_aspera' } 0..$#s;
    } else {
        my @s2 = split /;/, $s[$col];
        for (@s2) {
         print "$_\n"; # 输出结果用管道符传回perl
        }
    }
}
EOF

( [ ! -d $out ] ) && mkdir $out
parallel -j 2 ascp -P33001 -QT -l500m -i $openssh era-fasp@{} $out :::: tmp_dwn_list && rm tmp_dwn_list
