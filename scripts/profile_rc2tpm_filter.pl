#!/usr/bin/perl
# encoding: utf-8
# author : Jinxin Meng
# created date: 2023-10-09, 21:03:23
# modified date: 2023-11-19, 21:03:23
use warnings;
use strict;

die "Usage: perl $0 [rc_profile] [gene_length] [out_file] [tpm_threshold]\n  TPM of gene > tpm_threshold will be retained.\n" unless @ARGV eq 4;
my ($gene_rc, $gene_len, $out, $threshold) = @ARGV;
my (%len, @sum) = ();

open I, "<$gene_len";
while (<I>) {
    chomp;
    my @s = split/\s+/;
    $len{$s[0]} = $s[1];
}
close I;

open I, "<$gene_rc";
while (<I>) {
    chomp;
    next if /name/;
    my @s = split/\s+/;
    for my $i (1..$#s){
        $sum[$i] += $s[$i]/$len{$s[0]};
    }
}

open O, ">$out";
seek I, 0, 0;
while (<I>) {
    chomp;
    if (/name/) {
        print O "$_\n";
        next;
    }
    my @s = split/\s+/;
    my $f = 0;
    my @l = ();
    for my $i (1..$#s){
        $l[$i] = $s[$i]/$len{$s[0]}/$sum[$i]*1000000;
        $l[$i] = 0 if $l[$i] < $threshold;
        $f += $l[$i];
    }
    $l[0] = $s[0];
    print O join("\t",@l)."\n" unless $f eq 0;
}
close O;
