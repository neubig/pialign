#!/usr/bin/perl

use List::Util qw( max min sum );
use strict;
binmode STDOUT, ":utf8";

if(@ARGV < 3 or @ARGV > 5) {
    print "Usage: visualize.pl EFILE FFILE AFILE [ADDPOS] [REVERSE]\n";
    exit 1;
}

open EFILE, "<:utf8", $ARGV[0] or die "$ARGV[0]: $!";
open FFILE, "<:utf8", $ARGV[1] or die "$ARGV[1]: $!";
open AFILE, "<:utf8", $ARGV[2] or die "$ARGV[2]: $!";
my $EWIDTH = $ARGV[3] ? $ARGV[3] : 1;
my $FWIDTH = $ARGV[4] ? $ARGV[4] : 1;

my ($i, $j);
my (@actives, $actmax, $estr, $fstr, $astr);
while($estr = <EFILE> and $fstr = <FFILE> and $astr = <AFILE>) {
    chomp $estr; chomp $fstr; chomp $astr;
    # if(not $astr) { next; print "\n"; }
    my %active = map { my ($e,$f) = split(/-/); my $id = "".$e."-".$f; $id => 1 } split(/ /,$astr);
    my @e = split(/ /,$estr);
    my $elen = max( map { length } @e );
    my @f = split(/ /,$fstr);
    my $flen = max( map { length } @f );
    for(0 .. $flen) {
        print " " for(0 .. $elen*$EWIDTH);
        my $pos = $flen-$_;
        for(@f) {
            if($pos<length($_)) {
                print substr($_,length($_)-$pos-1,1);
            } else {
				print " " for( 1 .. $FWIDTH);
            }
        }
        print "\n";
    }
    foreach my $i (0 .. $#e) {
        print " " for(1 .. ($elen-length($e[$i]))*$EWIDTH);
        print "$e[$i] ";
        foreach my $j (0 .. $#f) {
            my $id = "$i-$j";
            print ($active{$id}?"X":".");
			print " " for(2 .. $FWIDTH);
        }
        print "\n";
    }
    print "\n";
}
