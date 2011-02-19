#!/usr/bin/perl

use List::Util qw( min max sum );

binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";
binmode STDERR, ":utf8";

if(@ARGV != 2) {
    print STDERR "Usage: lexicalize-pt.pl F2E_LEX E2F_LEX < unlex.pt > lex.pt\n";
    exit 1;
}

my $MIN = 1e-10;

# load the forward and backward lexical probabilities
open F2E, "<:utf8", $ARGV[0];
my %f2e = map { /^(.*)[ \t]([^ \t]*)$/; $1 => $2 } <F2E>; 
close F2E;
open E2F, "<:utf8", $ARGV[1];
my %e2f = map { /^(.*)[ \t]([^ \t]*)$/; $1 => $2 } <E2F>; 
close E2F;

sub lex {
    my ($e, $f, $s) = @_;
    my @es = split(/ /,$e);
    my @fs = split(/ /,$f);
    my $prob = 1;
    foreach my $ew (@es) {
        my $word = 0;
        foreach my $fw (@fs) {
            # print "s->{'$ew $fw'} = ".$s->{"$ew $fw"}."\n";
            $word += max($s->{"$ew $fw"}, $MIN);
        }
        $prob *= max($word/@fs,$MIN);
    }
    return $prob;
}

while(<STDIN>) {
    chomp;
    my ($f,$e,$s) = split(/ \|\|\| /);
    my @ss = split(/ /, $s); my $e1 = pop @ss;
    push @ss, lex($f,$e,\%e2f), lex($e,$f,\%f2e), $e1;
    print "$f ||| $e ||| @ss\n";
}
