#!/usr/bin/perl

use List::Util qw(min max sum);

if(@ARGV < 2 or @ARGV > 5) {
    print STDERR "Usage: aer.pl REF TEST [sent|total|points] ETEXT FTEXT\n";
    exit 1;
}
my $type = ($ARGV[2]?$ARGV[2]:"total");

open REF, "<:utf8", $ARGV[0] or die "$ARGV[0]: $!";
open TEST, "<:utf8", $ARGV[1] or die "$ARGV[1]: $!";
if($type eq "points") {
    $ARGV[3] and $ARGV[4] or die "Must specify a file when using points";
    open FTXT, "<:utf8", $ARGV[3] or die "$ARGV[3]: $!";
    open ETXT, "<:utf8", $ARGV[4] or die "$ARGV[4]: $!";
}

my ($tas, $tap, $ta, $ts);
while($ref = <REF> and $test = <TEST>) {
    chomp $ref; chomp $test;
    my ($cp,$cs,$cas,$cap,%rhash,$maxe,$maxf);
    for(split(/ /, $ref)) {
        /^([PS])-(\d+)-(\d+)$/ or die "bad ref alignment $_";
        $cs++ if $1 eq "S";
        $rhash{"$2-$3"} = $1;
        $maxe = max($maxe,$2);
        $maxf = max($maxf,$3);
    }
    my @tarr = split(/ /,$test);
    my @bad;
    for(@tarr) {
        my $val = $rhash{$_};
        if($val) {
            $cap++;
            $cas++ if($val eq 'S');
            delete $rhash{$_};
        } else {
            push @bad, "P-$_";
        }
    }
    while(my ($k,$v) = each(%rhash)) {
        push @bad, "R-$k" if $v eq "S";
    }
    my $crec = ($cas/$cs)*100;
    my $cprec = ($cap/(@tarr?@tarr:1))*100;
    my $cerr = (1-($cas+$cap)/($cs+@tarr))*100;
    printf "S=%i, EL=%i, FL=%i, R=%.2f%%, P=%.2f%%, E=%.2f%%\n", ++$sentid,$maxe,$maxf,$crec,$cprec, $cerr if $type eq "sent";
    $ts += $cs;
    $ta += @tarr;
    $tas += $cas;
    $tap += $cap;
    if($type eq "points") {
        my $esent = <ETXT>; chomp $esent; my @earr = split(/ /,$esent);
        my $fsent = <FTXT>; chomp $fsent; my @farr = split(/ /,$fsent);
        for(@bad) {
            my ($type,$f,$e) = split(/-/);
            print "$type\t$farr[$f]\t$earr[$e]\n";
        }
    }
}

my $crec = ($tas/$ts)*100;
my $cprec = ($tap/$ta)*100;
my $cerr = (1-($tas+$tap)/($ts+$ta))*100;
printf "Total: R=%.2f%%, P=%.2f%%, E=%.2f%%\n",$crec,$cprec, $cerr if $type ne "points";
