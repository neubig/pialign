#!/usr/bin/perl

use List::Util qw(sum);

open FFILE, "<:utf8", $ARGV[0] or die $!;
open EFILE, "<:utf8", $ARGV[1] or die $!;

my ($totf,$tote,$alf,$ale);
while($_ = <STDIN> and $es = <EFILE> and $fs = <FFILE>) {
    chomp; chomp $es; chomp $fs;
    my @esa = split(/ /,$es);
    my @fsa = split(/ /,$fs);
    my (@f,@e);
    for(split(/ /)) {
        my ($fv,$ev) = split(/-/);
        $f[$fv] = 1;
        $e[$ev] = 1;
    }
    $totf += @fsa;
    $alf += sum(@f);
    $tote += @esa;
    $ale += sum(@e);
}

printf("unaligned f=%.2f%%, e=%.2f%%\n", ($totf-$alf)/$totf*100, ($tote-$ale)/$tote*100);
