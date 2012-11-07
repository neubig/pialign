#!/usr/bin/perl

use strict;
use utf8;
use Getopt::Long;
use List::Util qw(sum min max shuffle);
binmode STDIN, ":utf8";
binmode STDOUT, ":utf8";
binmode STDERR, ":utf8";

my $TM_ORDER = 7;
my $TM_FEATS = 7;
my $TM_BINARY = 0;
my $LM_ORDER = 5;
my $DM_FILE = undef;
my $result = GetOptions (
    "tm-order=s" => \$TM_ORDER,   # The maximum length of a phrase in the TM
    "tm-binary=s" => \$TM_BINARY, # 1 if the TM is binary, otherwise 0
    "tm-feats=s" => \$TM_FEATS,   # 1 if the TM is binary, otherwise 0
    "dm-file=s" => \$DM_FILE,     # The distortion model file, if it exists
    "lm-order=s" => \$LM_ORDER,   # The n-gram order of the language model
); # flag

print STDERR "@ARGV\n";

if(@ARGV != 2) {
    print STDERR "Usage: $0 TM_FILE LM_FILE\n";
    exit 1;
}
my $TM_FILE = $ARGV[0];
my $LM_FILE = $ARGV[1];

# The weights of the distortion model
my $dm_weights = join("\n", map { 0.3 } (1 .. ($DM_FILE ? 7 : 1)));
my $tm_weights = join("\n", map { 0.2 } (1 .. $TM_ORDER));

print "
# input factors
[input-factors]
0

# mapping steps
[mapping]
0 T 0

# translation tables: table type (hierarchical(0), textual (0), binary (1)), source-factors, target-factors, number of scores, file 
[ttable-file]
$TM_BINARY 0 0 $TM_ORDER $TM_FILE

# no generation models, no generation-file section

# language models: type(srilm/irstlm), factors, order, file
[lmodel-file]
8 0 $LM_ORDER $LM_FILE

# limit on how many phrase translations e for each phrase f are loaded
# 0 = all elements loaded
[ttable-limit]
20
".
($DM_FILE ?
"# distortion (reordering) files
[distortion-file]
0-0 wbe-msd-bidirectional-fe-allff 6 $DM_FILE" :
"")
."

# distortion (reordering) weight
[weight-d]
$dm_weights

# language model weights
[weight-l]
0.5

# translation model weights
[weight-t]
$tm_weights

# word penalty
[weight-w]
-1
";
