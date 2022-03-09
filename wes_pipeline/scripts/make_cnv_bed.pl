#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';
use List::Util qw( min max );

#INPUT

my $infile=$ARGV[0];
my $outfile=$ARGV[1];
my @line=();

open B, "<$infile";
open C, ">$outfile";
while (<B>){
	chomp;
	@line = split;
	next if ($line[0] =~ m'#');
	if ($line[0] ne 'CHR') {
		my $chr='';
		my $start='';
		my $end='';
		my $type='';
		my @temp=();

		@temp = split ':', $line[2];
#		print @temp;
		$chr = $temp[2];
		$start = ((split '-', $temp[3])[0]);
		$end = ((split '-', $temp[3])[1]);
		if ($temp[1] eq 'GAIN') {
			$type="DUP";
		}
		elsif ($temp[1] eq 'LOSS') {
			$type="DEL";
		}
		print C "$chr\t$start\t$end\t$type\n"
	}
}

close B;
close C;