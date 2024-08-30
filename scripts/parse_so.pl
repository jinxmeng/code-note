#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-01-01, 12:07:31
# modified date: 2023-09-25, 12:07:31
use warnings;
use strict;
die "Usage: perl $0 [so_file_directory] [out_file]\n" unless @ARGV eq 2;
my ($in, $out) = @ARGV;
my @file = split /\s+/, readpipe("find $in -name '*so'");

open O, ">$out";
print O "name\tcontigs\ttotal_len\tgap\tmean_len\tN50\tN90\tmax_len\tmin_len\tGC\n";
for my $i (@file) {
    if ($i =~ /\/+(.*).so$/) {
        my $n = $1;
        open I, "<$i";
        my @line = ();
        while (<I>) {
            chomp;
            if ($_ =~ /\):\s+(\S+?)\s/) {
                push @line, $1;
            }
        }
        print O "$n\t".join("\t",@line)."\n";
    }
}
