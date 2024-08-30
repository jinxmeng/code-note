#!/usr/bin/perl
# encoding: utf-8
# author: Jinxin Meng
# created date: 2023-05-06, 12:01:05
# modified date: 2024-03-10, 17:50:35
use warnings;
use strict;
use POSIX;

die "Usage: perl $0 [CAZymes blast for genome, reg. exp. <sample>_<contigs>_<gene> ] [out_prefix]\n" unless @ARGV eq 2;
my ($in, $out) = @ARGV;
my ($x, %genome2CAZymes, %count, %all_CAZymes, %all_CAZyfams) = ("");

open I, "<$in" or die "Can't open $in: $!\n";
open O, ">${out}_gene2CAZymes"; 
while (<I>) {
    chomp;
    my @s = split /\s+/;
    if ($x ne $s[0]) {
        $_ =~ /(\S+?)(_\S+?)\s\S+?\|(\S+?)\s/;
        my ($genome, $gene, $CAZyme) = ($1, "$1$2", $3);
        print O "$genome\t$gene\t$CAZyme\n";
        push @{$genome2CAZymes{$genome}}, $CAZyme; 
    }
    $x = $s[0];    
}
close O;

for my $i (keys %genome2CAZymes) { # tidying CAZymes' names
    for my $j (@{$genome2CAZymes{$i}}) {
        my @s = split /\|/, $j;
        for my $k (@s) {
            $k =~ s/_\S+//g;
            $count{$i}{$k} += 1/($#s+1);
            $all_CAZymes{$k} += 1/($#s+1);
        }
    }
}

open O, ">${out}_CAZymes_statistics";
for my $i (sort keys %all_CAZymes) { # count all CAZymes
    my $n = int($all_CAZymes{$i}*1000)/1000;
    print O "$i\t$n\n";
    (my $x = $i) =~ s/\d+//g;
    $all_CAZyfams{$x} += $all_CAZymes{$i};
}
close O;

open O, ">${out}_CAZyfams_statistics";
for my $i (sort keys %all_CAZyfams) { # count all CAZyme families
    my $n = int($all_CAZyfams{$i}*1000)/1000;
    print O "$i\t$n\n";
}
close O;

my @genome = sort keys %count;
my @CAZyme = sort keys %all_CAZymes;
open O, ">${out}_genome2CAZymes_profile";
print O "name\t".join("\t", @genome)."\n";
for my $i (@CAZyme) {
    my @line = ();
    for my $j (@genome) {
        if (exists $count{$j}{$i}) {
            push @line, int($count{$j}{$i}*1000)/1000;
        } else {
            push @line, 0;
        }
    }
    print O "$i\t".join("\t", @line)."\n";
}
close O;
