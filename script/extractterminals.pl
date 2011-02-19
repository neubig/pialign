#!/usr/bin/perl

while(<STDIN>) {
    chomp;
    for(split(/ /)) {
        $vocab{$_}++ if not $_ =~ /^\[\]<>$/;
    }
}

for(sort { $vocab{$b} <=> $vocab{$a} } keys %vocab) {
    print "$_\t$vocab{$_}\n";
}
