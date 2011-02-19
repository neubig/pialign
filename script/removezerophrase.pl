#!/usr/bin/perl

use List::Util qw(min max);

my $max = ($ARGV[0]?$ARGV[0]:7);
my $min = 1;
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";

while(<STDIN>) {
    chomp;
    my @arr = split(/ \|\|\| /);
    die "Bad array @arr" if(@arr != 3);
    my @lens = map { scalar(split(/ /)) } @arr;
    print "$_\n" if min($lens[0],$lens[1])>=$min and max($lens[0],$lens[1])<=$max;
}
