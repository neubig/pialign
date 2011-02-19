#!/usr/bin/perl

use strict;
use List::Util qw(min max);
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";

my $REVERSE = ($ARGV[0] eq 'reverse');

my (%fe);

sub findmid {
    my $c = 0;
    for(my $i = 0; $i < @_; $i++) {
        if($_[$i] eq '(((') {
            while (($i < @_) and ($_[++$i] ne ')))')) { }
        }
        $c++ if($_[$i] =~ /^[\[<]$/);
        $c-- if($_[$i] =~ /^[>\]]$/);
        return $i+1 if not $c;
    }
    die "could not find middle of: @_";
}

my (%topl, %topr, %botl, %botr);
sub printphrases {
    my ($node, $totalf, $totale) = @_;
    # terminal
    my ($fs,$fl,$es,$el,$id);
    if($node->[0] == 0) {
        $id = $node->[3];
        ($fs,$fl,$es,$el) = ($node->[4],$node->[1],$node->[5],$node->[2]);
    }
    # forward
    else {
        my $left = printphrases($node->[3],$totalf,$totale);
        my $right = printphrases($node->[4],$totalf,$totale);
        my ($lf,$le) = split(/ \|\|\| /,$left);
        my ($rf,$re) = split(/ \|\|\| /,$right);
        my $newf = $lf.(($lf and $rf)?" ":"").$rf;
        my $newe = (($node->[0]==1)?$le:$re).(($le and $re)?" ":"").(($node->[0]==1)?$re:$le);
        $id = "$newf ||| $newe";
        ($fs,$fl,$es,$el) = ($node->[5],$node->[3]->[1]+$node->[4]->[1],$node->[6],$node->[3]->[2]+$node->[4]->[2]);
    }
    # print "printphrases(@$node, tf=$totalf, te=$totale, fs=$fs, fl=$fl, es=$es, el=$el, $id)\n";
    if($fl and $el) {
        my ($prev,$next) = ("other","other");
        # print "botr{".($fs-1) ."|".($es-1)."}\n";
        # print "botl{".($fs-1) ."|".($es+$el)."}\n";
        # print "topl{".($fs-$fl) ."|".($es+$el)."}\n";
        # print "topr{".($fs-$fl) ."|".($es-1)."}\n";
        if($botr{"".($fs-1) ."|".($es-1)}) {
            $prev = "mono";
        } elsif($botl{"".($fs-1) ."|".($es+$el)}) {
            $prev = "swap";
        }
        if($topl{"".($fs+$fl) ."|".($es+$el)}) {
            $next = "mono"; 
        } elsif($topr{"".($fs+$fl) ."|".($es-1)}) {
            $next = "swap";
        }
        my $str = $id;
        print "$str ||| $prev $next\n";
    }

    
    return $id;
}

# nodes
my @space;
sub markwords {
    my ($node, $fstart, $estart) = @_;
    push @$node, $fstart, $estart;
    # existing phrase
    # print join('',@space)."fstart=$fstart, estart=$estart, node1=".$node->[1].", node2=".$node->[2]."\n";
    push @space," ";
    if($node->[1] and $node->[2]) {
        my $fend = $fstart+$node->[1]-1;
        my $eend = $estart+$node->[2]-1;
        $topl{"$fstart|$estart"}++;
        $topr{"$fend|$estart"}++;
        $botl{"$fstart|$eend"}++;
        $botr{"$fend|$eend"}++;
    }
    # forward
    if($node->[0] == 1) {
        markwords($node->[3],$fstart,$estart);
        markwords($node->[4],$fstart+$node->[3]->[1],$estart+$node->[3]->[2]);
    }
    # backward
    elsif($node->[0] == -1) {
        markwords($node->[3],$fstart,$estart+$node->[4]->[2]);
        markwords($node->[4],$fstart+$node->[3]->[1],$estart);
    }
    pop @space;
}

# returns NT: [ type, f length, e length, left node, right node, f start, e start ]
#      or  T: [ type, f length, e length, str, f start, e start]
#  type: 0=terminal, 1=forward, -1=backwards
sub buildlens {
    return 0 if(!@_);
    if($_[0] eq '(((') {
        $_[-1] eq ')))' or die "bad string @_";
        my ($f,$e) = split(/ \|\|\| /,join(' ',@_[1 .. $#_-1]));
        my $flen = scalar(split(/ /,$f));
        my $elen = scalar(split(/ /,$e));
        return [ 0, $flen, $elen, "$f ||| $e" ];
    }
    my ($s, @m) = @_;
    my $e = pop @m;
    (($s eq '[') and ($e eq ']')) or (($s eq '<') and ($e eq '>')) or die "bad arr (s=$s, e=$e): @_";
    my $mid = findmid(@m);
    my $left = buildlens(@m[0 .. $mid-1]);
    my $right = buildlens(@m[$mid .. $#m]);
    return [ (($s eq '[')?1:-1), $left->[1] + $right->[1], $left->[2] + $right->[2], $left, $right ];
}

while(<STDIN>) {
    chomp;
    s/{ //g; s/ }//g;
    my $root = buildlens(split(/ /));
    if($root) {
        %botr = ("-1|-1" => 1);
        %botl = ("-1|".$root->[2] => 1);
        %topl = ("".$root->[1]."|".$root->[2] => 1);
        %topr = ("".$root->[1]."|-1" => 1);
        markwords($root,0,0);
        # print "topl=".join(' ',keys %topl)."\n";
        # print "botl=".join(' ',keys %botl)."\n";
        # print "topr=".join(' ',keys %topr)."\n";
        # print "botr=".join(' ',keys %botr)."\n";
        printphrases($root,$root->[1],$root->[2]);
    }
}
