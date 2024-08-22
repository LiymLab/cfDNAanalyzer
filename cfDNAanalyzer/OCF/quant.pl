#!/usr/bin/perl

#
# Author: Ahfyth
# Date  : Dec 12, 2108
#
# This program is designed to quantify the orientation-aware cfDNA ends
# in tissue-specific open chromatin regions (OCF values)
#
# OCF: Orientation-aware CfDNA Fragmentation
#

use strict;
use warnings;

if( $#ARGV < 0 )
{
	print STDERR "\nUsage: $0 <in.spc.sync.end> [peak=60]\n\n";
	exit 2;
}

## Calculate the OCF (Orientation-aware CfDNA Fragmentation) value
my $peak = $ARGV[1] || 60;
my $bin  = 10;

## Definition 1
my $trueends   = 0;	## neg (-peak-bin to -peak+bin) + pos (peak-bin to peak+bin)
my $background = 0;	## neg (peak-bin to peak+bin) + pos (-peak-bin to -peak+bin)

open IN, "$ARGV[0]" or die( "$!" );
<IN>;	#head
###Locus	#left	%%Left	#Right	%%Right
while( <IN> )
{
	chomp;
	my @l = split /\t/;
	if( $l[0]>=-$peak-$bin && $l[0]<=-$peak+$bin )
	{
		$trueends   += $l[4];
		$background += $l[2];
	}
	elsif( $l[0]>=$peak-$bin && $l[0]<=$peak+$bin )
	{
		$trueends   += $l[2];
		$background += $l[4];
	}

	last if $l[0] > $peak+$bin;
}
close IN;

print join("\t", $ARGV[0], $trueends-$background), "\n";

