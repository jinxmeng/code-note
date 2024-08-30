#!/usr/bin/perl
# encoding: utf-8
# author : Jinxin Meng
# created date: 2023-10-09, 21:06:53
# modified date: 2023-11-19, 21:06:53
use warnings;
use strict;

die "Usage: perl $0 [gene_to_KO] [gene_profile] [out_file] [threshold]\n  TPM of KO > tpm_threshold will be retained.\n" unless @ARGV eq 4;
my ($gene2KO, $gene_tpm, $out, $threshold) = @ARGV;
my (%KO, %sum) = ();

open I, "<$gene2KO";
while (<I>) {
    chomp;
    my @s = split /\s+/;
    $KO{$s[0]} = $s[1];
} 

open I, "<$gene_tpm";
open O, ">$out";
my $line = <I>;
print O "$line";
while(<I>){
    chomp;
    my @s = split /\t/;
    if (exists $KO{$s[0]}) {
        for my $i (1..$#s) {
            @{$sum{$KO{$s[0]}}}[$i-1] += $s[$i];
        }
    }
}
close I;

for my $i (sort keys %sum) {
    my $flag = 0;
    for my $j (0..$#{$sum{$i}}) {
        @{$sum{$i}}[$j] = 0 if @{$sum{$i}}[$j] < $threshold;
        $flag += @{$sum{$i}}[$j];
    }
    print O "$i\t".join("\t",@{$sum{$i}})."\n" unless $flag eq 0;
} 
close O;
