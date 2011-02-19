#!/usr/bin/perl

use List::Util qw(min max sum);
use strict;

sub quasi_consecutive {
    my $small = shift;
    my $large = shift;
    my %tp = %{shift(@_)};
    my @f2e = @{shift(@_)};
    for($small .. $large) {
        return 0 if ($f2e[$_] and not exists $tp{$_});
    }
    return 1;
}

while(<STDIN>) {
    chomp;
    s/[SP]-//g;
    my ($I,$J,@e2f,@f2e);
    # save in appropriate format
    for(split(/ /)) {
        my ($e,$f) = split(/-/);
        # print "e-f= $e-$f\n";
        $I = max($I,$e);
        $J = max($J,$f);
        $e2f[$e] = [] if not exists $e2f[$e];
        push @{$e2f[$e]}, $f;
        $f2e[$f] = [] if not exists $f2e[$f];
        push @{$f2e[$f]}, $e;
    }
    
    # phrase-extract
    my @out;
    foreach my $i1 (0 .. $I) {
        my %tp;
        foreach my $i2 ($i1 .. $I) {
            if($e2f[$i2]) { for (@{$e2f[$i2]}) { $tp{$_}++; } }
            my $j1 = min(keys %tp);
            my $j2 = max(keys %tp);
            if(quasi_consecutive($j1, $j2, \%tp, \@f2e)) {
                my %sp;
                foreach my $j ($j1 .. $j2) {
                    if($f2e[$j]) { for (@{$f2e[$j]}) { $sp{$_}++; } }
                }
                if(min(keys %sp) >= $i1 and max(keys %sp) <= $i2) {
                    while($j1 >= 0) {
                        my $jp = $j2;
                        while($jp <= $J) {
                            push @out, "$i1-".($i2+1)."-$j1-".($jp+1);
                            $jp++;
                            last if $f2e[$jp];
                        }
                        $j1--;
                        last if $f2e[$j1];
                    }
                }
            }
        }
    }
    print "@out\n";
}
