#!/usr/bin/perl

#
# Author: Ahfyth
# Date  : Dec 12, 2018
#
# This program is designed to extract the orientation-aware cfDNA ends information
# around tissue-specific open chromatin regions
#

use strict;
use warnings;

if( $#ARGV < 1 )
{
	print STDERR "\nUsage: $0 <in.ol> <out.prefix> [extended=1000]\n";
	exit 2;
}

my $ext = $ARGV[2] || 1000;

my (%left, %right);
my %total;
open IN, "$ARGV[0]" or die( "$!" );
while( <IN> )
{
	chomp;
	my @l = split /\t/;	##chr1	1014459	1016459	LABEL	chr1	1015296	1015412	116
	next if $l[4] eq '.' || $l[0]!~/^chr\d/;
	if( $l[5] >= $l[1] )
	{
		my $s = $l[5] - $l[1];
		$left{ $l[3] }->[ $s ] ++;
		$total{ $l[3] }->{S} ++;
	}
	if( $l[6] <= $l[2] )
	{
		my $e = $l[6]-($l[1]+1);
		$right{ $l[3] }->[ $e ] ++;
		$total{ $l[3] }->{E} ++;
	}
}
close IN;

# 2024.8.6 zjp modified
#foreach my $k ( keys %total )
#{
#	open OUT, ">$ARGV[1].$k.sync.end" or die( "$!" );

#	my $le = $left{$k};
#	my $re = $right{$k};
#	my $ts = $total{$k}->{S}/10000;
#	my $te = $total{$k}->{E}/10000;
#	my $num = $ext*2-1;

#	print OUT "#Locus\t#left\t%%Left\t#Right\t%%Right\n";
#	foreach my $k ( 0..$num )
#	{
#		my $l = $le->[$k] || 0;
#		my $r = $re->[$k] || 0;
#		print OUT join("\t", $k-$ext, $l, $l/$ts, $r, $r/$te), "\n";
#	}
#	close OUT;
#}

foreach my $k ( keys %total )
{
    my $ts = $total{$k}->{S} // 0;  # 使用 // 确保变量有默认值
    my $te = $total{$k}->{E} // 0;

    # 防止除数为 0
    $ts = 1 if $ts == 0;
    $te = 1 if $te == 0;

    $ts /= 10000;
    $te /= 10000;

    open OUT, ">$ARGV[1].$k.sync.end" or die( "$!" );

    my $le = $left{$k};
    my $re = $right{$k};
    my $num = $ext*2-1;

    print OUT "#Locus\t#left\t%%Left\t#Right\t%%Right\n";
    foreach my $k ( 0..$num )
    {
        my $l = $le->[$k] || 0;
        my $r = $re->[$k] || 0;
        print OUT join("\t", $k-$ext, $l, $l/$ts, $r, $r/$te), "\n";
    }
    close OUT;
}

# end 