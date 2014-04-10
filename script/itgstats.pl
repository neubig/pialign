#!/usr/bin/perl

use strict;
use List::Util qw(min max);
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";

# a script to calculate a number of statistics about ITGs
# ./itgstats.pl TYPE [MAX_LENGTH=7] [OPTIONS]
#
# TYPE is the type of processing to perform
#
# The types of processing that can be performed are:
#  merge: Take an ITG derivation file and merge into phrase pairs
#  block: Take an ITG derivation file and merge into 1-many blocks
#  lex:   Calculate the lexical reordering probabilities of each phrase
#         pair. OPTION is the smoothing pseudo-count to add to each
#         type of reordering (default 0.5)
#  phrase: Build a phrase table from the alignments using maximum
#          likelihood estimation
#  align/palign/balign: Output word alignments in GIZA++ format using
#         1-1 words, 1-many blocks, or many-many phrases respectively

if(@ARGV < 1) {
    print STDERR "Usage: itgstats.pl merge|block|lex|phrase|align|palign|balign|talign|salign [MAX_LENGTH=7] [OPTIONS]\n";
    exit 1;
}

my $TYPE = $ARGV[0];
my $MAXLEN = $ARGV[1] ? $ARGV[1] : 7;
my $OPT1 = $ARGV[2] ? $ARGV[2] : 0.5;
my $MINLEN = 1;
my $MERGE = (($TYPE eq "merge") or ($TYPE eq "palign"));
my $BLOCK = (($TYPE eq "block") or ($TYPE eq "balign"));
my $TRGONE = $TYPE eq "talign";
my $SRCONE = $TYPE eq "salign";

my (%fe, %phrases);

sub findmid {
    my $c = 0;
    for(my $i = 0; $i < @_; $i++) {
        if($_[$i] eq '(((') {
            while (($i < @_) and ($_[++$i] ne ')))')) { }
        }
        $c++ if($_[$i] =~ /^[\[<{]$/);
        $c-- if($_[$i] =~ /^[>\]}]$/);
        return $i+1 if not $c;
    }
    die "could not find middle of: @_";
}


# print the tree for debugging
my @spaces;
sub printtree {
    my $tree = shift;
    return if not $tree;
    push @spaces, " ";
    print join("",@spaces).join(',',@$tree)."\n";
    if($tree->[6]) {
        printtree($tree->[6]); printtree($tree->[7]);
    }
    pop @spaces;
}

# find the starting points given the tree
sub buildstarts {
    my ($tree, $fstart, $estart) = @_;
    return if not $tree;
    $tree->[1] = $fstart; $tree->[2] = $estart;
    if($tree->[0] == 1) {
        buildstarts($tree->[6],$fstart,$estart);
        buildstarts($tree->[7],$fstart+$tree->[6]->[3],$estart+$tree->[6]->[4]);
    } elsif($tree->[0] == -1) {
        buildstarts($tree->[6],$fstart,$estart+$tree->[7]->[4]);
        buildstarts($tree->[7],$fstart+$tree->[6]->[3],$estart);
    }
}

# returns: [ type, fstrt, estrt, flen, elen, phrase, lnode, rnode]
#  type: 0=terminal, 1=forward, -1=backwards
sub buildtree {
    my $inphrase = shift;
    return 0 if(!@_);
    # terminals
    if($_[0] eq '(((') {
        $_[-1] eq ')))' or die "bad string @_";
        my ($f,$e) = split(/ \|\|\| /,join(' ',@_[1 .. $#_-1]));
        my $flen = scalar(split(/ /,$f));
        my $elen = scalar(split(/ /,$e));
        return [ 0, -1, -1, $flen, $elen, "$f ||| $e" ];
    }
    # if a phrase, set inphrase to be true
    if($_[0] eq '{') {
        $_[-1] eq '}' or die "bad string @_";
        $inphrase = 1;
        @_ = @_[1 .. $#_-1];
    }
    # do forward or backward
    my ($s, @m) = @_;
    my $e = pop @m;
    (($s eq '[') and ($e eq ']')) or (($s eq '<') and ($e eq '>')) or die "bad arr (s=$s, e=$e): @_";
    my $mid = findmid(@m);
    my $l = buildtree($inphrase, @m[0 .. $mid-1]);
    my $r = buildtree($inphrase, @m[$mid .. $#m]);
    my ($lf,$le) = split(/ \|\|\| /,$l->[5]);
    my ($rf,$re) = split(/ \|\|\| /,$r->[5]);
    my $t = (($s eq '[')?1:-1);
    my $phrase = $lf.((length($lf) and length($rf))?" ":"").$rf." ||| ".($t==1?$le:$re).((length($le) and length($re))?" ":"").($t==1?$re:$le);
    my $ret =  [ $t, -1, -1, $l->[3]+$r->[3], $l->[4]+$r->[4], $phrase, $l, $r ];
    # throw away children if this is to be merged
    if($inphrase and ($MERGE or ($BLOCK and (min($ret->[3],$ret->[4]) <= 1)) or ($TRGONE and ($ret->[4] <= 1)) or ($SRCONE and ($ret->[3] <= 1)) or max($ret->[3],$ret->[4]) == 1)) {
        pop @$ret; pop @$ret; 
        $ret->[0] = 0;
    }
    return $ret;
}

# mark the corners for lexical probabilities
sub markcorners {
    my ($node, $corners) = @_;
    if($node->[3]*$node->[4]) {
        $corners->{"tl".$node->[1]."|".$node->[2]}++;
        $corners->{"bl".$node->[1]."|".($node->[2]+$node->[4]-1)}++;
        $corners->{"tr".($node->[1]+$node->[3]-1)."|".$node->[2]}++;
        $corners->{"br".($node->[1]+$node->[3]-1)."|".($node->[2]+$node->[4]-1)}++;
    }
    if($node->[0]) {
        markcorners($node->[6],$corners);
        markcorners($node->[7],$corners);
    }
}

# count the lexical reordering frequencies
sub countlex {
    my ($node, $corners) = @_;
    if($node->[0]) {
        countlex($node->[6],$corners);
        countlex($node->[7],$corners);
    }
    if(min($node->[3],$node->[4])>=$MINLEN and max($node->[3],$node->[4])<=$MAXLEN) {
        if(not $phrases{$node->[5]}) {
            my @arr = map { $OPT1 } (1 .. 6);
            $phrases{$node->[5]} = \@arr;
        }
        my $arr = $phrases{$node->[5]};
        if($corners->{"br".($node->[1]-1)."|".($node->[2]-1)})             { $arr->[0]++; }
        elsif($corners->{"tr".($node->[1]-1)."|".($node->[2]+$node->[4])}) { $arr->[1]++; }
        else                                                               { $arr->[2]++; }
        if($corners->{"tl".($node->[1]+$node->[3])."|".($node->[2]+$node->[4])}) { $arr->[3]++; }
        elsif($corners->{"bl".($node->[1]+$node->[3])."|".($node->[2]-1)})       { $arr->[4]++; }
        else                                                                     { $arr->[5]++; }
    }
}


# print the terminal alignments
my $first = 1;
sub printalign {
    my ($node) = @_;
    return if not $node;
    if($node->[0]) {
        printalign($node->[6]);
        printalign($node->[7]);
    } else {
        for(my $i = 0; $i < $node->[3]; $i++) {
            for(my $j = 0; $j < $node->[4]; $j++) {
                if(not $first) {
                    print " ";
                }
                $first = 0;
                print "".($i+$node->[1])."-".($j+$node->[2])."";
            }
        }
    }
}

# print the ITG in the same format as the input
sub printitg {
    my $node = shift;
    return if(not $node);
    if($node->[0]) {
        print ($node->[0]==1?"[ ":"< ");
        printitg($node->[6]);
        print " ";
        printitg($node->[7]);
        print ($node->[0]==1?" ]":" >");
    } else {
        print "((( ".$node->[5]." )))";
    }
}

sub countphrase {
    my $node = shift;
    if($node) {
        $phrases{$node->[5]}++;
        if($node->[0]) {
            countphrase($node->[6]);
            countphrase($node->[7]);
        }
    }
}

while(<STDIN>) {
    chomp;
    my $root = buildtree(0, split(/ /));
    buildstarts($root,0,0);
    # printtree($root);
    if($TYPE eq "lex") {
        if($root) {
            my %corners = ( "br-1|-1" => 1, "bl".$root->[3]."|-1" => 1,
                    "tr-1|".$root->[4] => 1, "tl".$root->[3]."|".$root->[4] => 1 ); 
            markcorners($root,\%corners);
            # print join(' ',sort keys(%corners))."\n";
            countlex($root,\%corners);
        }
    } elsif($TYPE =~ /(merge|block)/) {
        printitg($root);
        print "\n";
    } elsif($TYPE eq "phrase") {
        countphrase($root);
    } elsif($TYPE =~ /^[pbts]?align$/) {
        $first = 1;
        printalign($root);
        print "\n";
    } else {
        die "Invalid processing type $TYPE";
    }
}

if($TYPE eq "phrase") {
    my (%es, %fs, $tot);
    while(my ($k,$v) = each(%phrases)) {
        my($f,$e) = split(/ \|\|\| /,$k);
        $es{$e} += $v;
        $fs{$f} += $v;
        $tot += $v;
    }
    while(my ($k,$v) = each(%phrases)) {
        my($f,$e) = split(/ \|\|\| /,$k);
        my $flen = split(/ /,$f);
        my $elen = split(/ /,$e);
        if(min($flen,$elen)>=$MINLEN and max($flen,$elen) <= $MAXLEN) {
            print "$k ||| ".($v/$es{$e})." ".($v/$fs{$f})." ".($v/$tot)." 1 2.718\n";
        }
    }
}

if($TYPE eq "lex") {
    while(my ($k,$v) = each(%phrases)) {
        my @s = ($v->[0]+$v->[1]+$v->[2],$v->[3]+$v->[4]+$v->[5]);
        print "$k |||";
        for(0 .. 5) { print " ".($v->[$_]/$s[int($_/3)]); };
        print "\n";
    }
}
