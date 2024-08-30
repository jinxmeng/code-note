#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-10-02, 01:49:47
# modified date: 2024-04-28, 23:39:57

use warnings;
use strict;
use Getopt::Long;

my $usage;
$usage .= "Usage: perl $0 [OPTIONS] [Cdb.csv] [out_file]\n";
$usage .= "  OPTIONS:\n";
$usage .= "    -ckm  : based on quality score provided by checkm\n";
$usage .= "    -ckm2 : based on quality score provided by checkm2\n";
$usage .= "    -len  : based on genome len, tab-delimited file with field name and len\n";
$usage .= "  EXAMPLE:\n";
$usage .= "    perl parse_dRep.pl Cdb.csv genome_cls_res\n";
$usage .= "    perl parse_dRep.pl -ckm2 quality_report.tsv Cdb.csv genome_cls_res\n";

my ($ckm, $ckm2, $len) = ("", "", "");
&GetOptions(
    "ckm=s" => \$ckm,
    "ckm2=s" => \$ckm2,
    "len=s" => \$len);

die "$usage" unless @ARGV eq 2;
my ($cdb, $out) = @ARGV;

open C, "<$cdb" or die "Can't open $cdb: $!\n";
open O, ">$out" or die "Can't open $out: $!\n";

my (%cls, @name) = (); 
while (<C>) {
    chomp;
    next if (/primary_cluster/);
    my @s = split /,/;
    $s[0] =~ s/.(fa|fna)//;
    push @{$cls{$s[1]}}, $s[0];
    push @name, $s[0];
}

if ($ckm ne "") {
    my %qs;
    open I, "$ckm" or die "Can't open $ckm: $!\n";
    while (<I>) {
        chomp;
        next if /Completeness/;
        my @s = split /\t/;
        $qs{$s[0]} = $s[11] - 5 * $s[12] if (grep {$_ eq $s[0]} @name);
    }
    select_rep(\%cls, \%qs); # 调用子程序并传递哈希的引用
} elsif ($ckm2 ne "") {
    my %qs;
    open I,"<$ckm2" or die "Can't open $ckm2: $!\n";
    while (<I>) {
        chomp;
        next if /Completeness/;
        my @s = split /\t/;
        $qs{$s[0]} = $s[1] - 5 * $s[2] if (grep {$_ eq $s[0]} @name);
    }
    select_rep(\%cls, \%qs); # 调用子程序并传递哈希的引用
} elsif ($len ne "") {
    my %glen;
    open I, "<$len" or die "Can't open $len: $!\n";
    while (<I>) {
        chomp;
        my @s = split/\t/;
        $glen{$s[0]} = $s[1];
    }
    select_rep(\%cls, \%glen); # 调用子程序并传递哈希的引用
} else {
    for my $i (sort keys %cls) {
        my $count = scalar @{$cls{$i}};
        print O "cluster.$i\t$count\t".join(",", @{$cls{$i}})."\n";     
    }
}

close I;
close O;

sub select_rep {
    my ($cls, $x) = @_; # 获取哈希的引用
    my %cls = %{$cls}; # 将哈希引用解引用为哈希
    my %x = %{$x}; # 将哈希引用解引用为哈希
    for my $i (sort keys %cls) {
        my $count = scalar @{$cls{$i}};
        if ($count eq 1) {
            print O "cluster.$i\t1\t@{$cls{$i}}\t@{$cls{$i}}\n";
        } else {
            my $max = -10000;
            my $rep = "";
            for my $j (@{$cls{$i}}) {
                if ($x{$j} >= $max) {
                    $max = $x{$j};
                    $rep = $j;
                }
            }
            print O "cluster.$i\t$count\t$rep\t".join(",", @{$cls{$i}})."\n";
        }   
    }
}
