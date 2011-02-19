#!/usr/bin/perl

while(<STDIN>) {
    chomp;
    my @arr = split(/ /);
    for(0 .. $#arr) {
        my ($e,$f) = split(/-/,$arr[$_]);
        $arr[$_] = "$f-$e";
    }
    print "@arr\n";
}
