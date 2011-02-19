#!/usr/bin/perl

my @likelihoods;
foreach my $f (0 .. $#ARGV) {
    print ",$ARGV[$f]";
    my $iter = 0;
    open FILE, "<:utf8", $ARGV[$f] or die $!;
    while(<FILE>) {
        chomp;
        if(/Likelihood=([^ ,].*)/) {
            $likelihoods[$iter] = [] if not $likelihoods[$iter];
            $likelihoods[$iter]->[$f] = $1;
            $iter++;
        }
    }
    close FILE;
}

print "\n";
for(0 .. $#likelihoods) {
    print "". $_+1 .",".join(",", @{$likelihoods[$_]})."\n";
}
