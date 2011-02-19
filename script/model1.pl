#!/usr/bin/perl

use strict;
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";

my $NUM_ITERS = 20;
my $TERMINATE = 1/1000;
my $PACK = "LL";
my $CUTOFF = 1e-7;

if(@ARGV != 2) {
    print "Usage: model1.pl FFILE EFILE\n";
    exit 1;
}

# get an id for a word, add it if it doesn't exists
#  getid(word,id_hash,id_vector)
sub getid {
    my $s = shift;
    my $h = shift;
    my $i = shift;
    if(not exists $h->{$s}) {
        $h->{$s} = @$i;
        push @$i, $s;
    }
    return $h->{$s};
}

# load a corpus
#  loadcorp(filename)
#  returns [array,hashmap]
sub loadcorp {
    my (@arr, %ids, @syms);
    &getid("NULL", \%ids, \@syms);
    open FILE, "<:utf8", $_[0] or die $!;
    while(<FILE>) {
        chomp;
        my @sent = map { &getid($_, \%ids, \@syms) } split(/ /);
        push @arr, \@sent;
    }
    close FILE;
    return ( \@arr, \@syms );
}

my ( $fcorp, $fsyms ) = &loadcorp($ARGV[0]);
my ( $ecorp, $esyms ) = &loadcorp($ARGV[1]);

print STDERR "loaded ".scalar(@$fcorp)." sentences\n";

die "Corpus files are of different sizes" if(@$fcorp != @$ecorp);

# add null word
unshift @$_, 0 for(@$ecorp);

# initialize to uniform
my $fsize = @$fsyms-1;
my $esize = @$esyms-1;
my (%t,$f,$e);
for(0 .. @$ecorp-1) {
    my @earr = @{$ecorp->[$_]};
    my @farr = @{$fcorp->[$_]};
    foreach $f (@farr) {
        foreach $e (@earr) {
            $t{pack($PACK,$f,$e)} = 1/$fsize;
        }
    }
}

# train t
my $lastl;
foreach my $iter ( 1 .. $NUM_ITERS ) {
    my (%count, @total, $l);
    # E step
    for(0 .. @$ecorp-1) {
        my @earr = @{$ecorp->[$_]};
        my @farr = @{$fcorp->[$_]};
        my @stotal;
        foreach $f (@farr) {
            foreach $e (@earr) {
                my $id = pack($PACK,$f,$e);
                $stotal[$f] += $t{$id};
            }
            $l += log($stotal[$f]/@earr);
        }
        foreach $f (@farr) {
            foreach $e (@earr) {
                my $id = pack($PACK,$f,$e);
                $count{$id} += $t{$id}/$stotal[$f];
                $total[$e] += $t{$id}/$stotal[$f];
            }
        }
    }
    # M step
    foreach my $id ( keys %t ) {
        my ($f,$e) = unpack($PACK, $id);
        $t{$id} = $count{$id}/$total[$e];
    }
    print STDERR "Iter $iter: l=$l\n";
    last if($lastl and ($lastl-$l)>$l*$TERMINATE);
    $lastl = $l;
}

for(sort keys %t) {
    if($t{$_} > $CUTOFF) {
        ($f,$e) = unpack($PACK, $_);
        printf("%s %s %0.7g\n", $fsyms->[$f], $esyms->[$e], $t{$_});
    }
}
