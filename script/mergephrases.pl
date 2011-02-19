#!/usr/bin/perl

use strict;
binmode STDOUT, ":utf8";

my (%fes,@fs,@es);
my (%ids, @strs);

sub getid {
    $_ = shift;
    if(not exists $ids{$_}) {
        $ids{$_} = scalar(keys %ids);
        push @strs,$_;
    }
    return $ids{$_};
}

for(@ARGV) {
    print STDERR "Processing $_\n";
    open FILE, "<:utf8", $_ or die $!;
    while(<FILE>) {
        chomp;
        my @arr = split(/ \|\|\| /);
        @arr == 3 or die "bad line $_\n";
        my @scores = split(/ /, $arr[2]);
        my $f = getid($arr[0]);
        my $e = getid($arr[1]);
        my $fe = "$f-$e";
        $fes{$fe} = [0,0] if not $fes{$fe};
        $fes{$fe}->[0] += $scores[2];
        $fes{$fe}->[1] += $scores[3];
        $fs[$f] += $scores[2];
        $es[$e] += $scores[2];
    }
    close FILE;
}

my $exp1 = exp(1);
while(my ($k,$v) = each(%fes)) {
    my ($f,$e) = split(/-/,$k);
    printf ("%s ||| %s ||| %.6g %.6g %.6g %.6g %.6g\n", $strs[$f], $strs[$e], ($v->[0]/$es[$e]), ($v->[0]/$fs[$f]), ($v->[0]/@ARGV), ($v->[1]/@ARGV), $exp1);
}
