#!/usr/bin/perl
# encoding: utf-8
# author : Jinxin Meng
# created date: 2023-10-10, 21:05:16
# modified date: 2023-12-07, 22:20:01
use warnings;
use strict;

die "Usage: perl $0 [gene_to_CAZyme (AA5|AA3_1 ..)] [gene_tpm] [out_file] [threshold]\n  TPM of CAZyme > tpm_threshold will be retained.\n" unless @ARGV eq 4;
my ($gene2enz, $gene_tpm, $out, $threshold) = @ARGV;
my (@name, %enz, %perc, %sum) = ();

open I, "<$gene2enz";
while (<I>) {
    chomp;
    my @s = split /\s+/;
    my @s2 = split /\|/, $s[1];
    # for my $i (@s2) {
    #     $i =~ s/_\d+//g;
    # }
    @{$enz{$s[0]}} = @s2;
    $perc{$s[0]} = scalar @s2;
    for my $i (@s2) {
        if (! grep { $_ eq $i } @name) {
            push @name, $i;
        }
    }
}

open I, "<$gene_tpm";
open O, ">$out"; 
my $h = <I>; 
print O "$h";
while(<I>){
    chomp;
    my @s = split /\s+/;
    if (exists $enz{$s[0]}) {
        for my $i (@{$enz{$s[0]}}) {
            for my $j (1..$#s) {
                @{$sum{$i}}[$j-1] += $s[$j] / $perc{$s[0]};
            }
        }
    }
}

for my $i (sort keys %sum) {
    my $flag = 0;
    for my $j (0..$#{$sum{$i}}) {
        @{$sum{$i}}[$j] = 0 if @{$sum{$i}}[$j] < $threshold;
        $flag += @{$sum{$i}}[$j];
    }
    print O "$i\t".join("\t",@{$sum{$i}})."\n" unless $flag eq 0;
} 
close O;
