#!/usr/bin/perl

my %nums = ("r"=>0, "i"=>1,"b"=>2,"g"=>3,"f"=>4,"s"=>5);
my $idx=1;
print "Iter, Total, Rem, Init, Base, Gen, For, Samp, Rem%, Init%, Base%, Gen%, For%, Samp%\n";
while(<STDIN>) {
    chomp;
    if(/Time=([0-9\.]*)s \((.*)\)/) {
        my @arr;
        $arr[0] = $1;
        $_ = $2;
        while(/([a-z])=([0-9\.]*)/g) {
            $arr[$nums{$1}+1] = $2;
            $arr[$nums{$1}+7] = (100*$2/$arr[0])."%";
        }
        print "$idx,".join(',',@arr)."\n";
        $idx++;
    }
}
