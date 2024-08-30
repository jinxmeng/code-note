#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-05-06, 12:01:05
# modified date: 2024-03-10, 21:16:21
use warnings;
use strict;
use POSIX;

die "Usage: perl $0 [CAZymes blast for geneset] [out_prefix]\n" unless @ARGV eq 2;
my ($in, $out) = @ARGV;
my ($x, %all_CAZymes, %all_CAZyfams) = ("");

open I, "<$in" or die "Can't open $in: $!\n";
open O, ">${out}_gene2CAZymes"; 
while (<I>) {
    chomp;
    my @s = split /\s+/;
    if ($x ne $s[0]) {
        $_ =~ /(\S+?)\s\S+?\|(\S+?)\s/;
        my ($gene, $CAZyme) = ($1, $2);
        print O "$gene\t$CAZyme\n";
        my @s = split/\|/, $CAZyme;
        for my $i (@s) {
            $i =~ s/_\S+//g;
            $all_CAZymes{$i} += 1/($#s+1);
        }
    }
    $x = $s[0];    
}
close O;

open O, ">${out}_CAZymes_statistics";
for my $i (sort keys %all_CAZymes) { # count all CAZymes
    my $n = int($all_CAZymes{$i}*1000)/1000;
    print O "$i\t$n\n";
    if ($i =~ /\S+\.\S+\.\S+/) {
        $all_CAZyfams{"Other"} += $all_CAZymes{$i};
    } else {
        (my $x = $i) =~ s/\d+//g;
        $all_CAZyfams{$x} += $all_CAZymes{$i};
    }
}
close O;

open O, ">${out}_CAZyfams_statistics";
for my $i (sort keys %all_CAZyfams) { # count all CAZyme families
    my $n = int($all_CAZyfams{$i}*1000)/1000;
    print O "$i\t$n\n";
}
close O;
