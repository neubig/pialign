#!/usr/bin/perl

my ($curr,@arr);
while(<STDIN>) {
    chomp;
    my ($sid,$e,$f,$t) = split(/ /);
    if($sid ne $curr) {
        print "@arr\n" if(@arr);
        @arr = ();
        $curr = $sid;
    }
    push @arr, "$t-".($e-1)."-".($f-1);
}
print "@arr\n" if(@arr);
