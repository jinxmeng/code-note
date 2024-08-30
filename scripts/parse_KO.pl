#!/usr/bin/perl
use strict;

die "Usage: perl $0 <ko00001.keg> <STDOUT>\n" unless @ARGV == 1;
my $in_f = @ARGV[0];

open IN,"<$in_f";

my ($a,$b,$c,$d) = ();
while(<IN>) {
    chomp;

    if(/^A/) {

        $a = join("\t",split(/\s/,$_,2))

    } elsif (/^B\s+(\d+)\s+(.*)/) {

        $b = join("\t","B".$1, $2);

    } elsif (/^C\s+(\d+)\s+(.*)/) {

        $c = join("\t","C".$1, $2);

    } elsif (/^D\s+K(\d+)\s+(.*)/) {

        $d = join("\t","K".$1, $2);
        print STDOUT "$a\t$b\t$c\t$d\n";

    } else { next }
}

close(IN);
############################################################
##print STDERR "Program End...\n";gg=G
