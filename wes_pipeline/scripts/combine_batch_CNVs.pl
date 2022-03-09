#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $infile = $ARGV[0];
my $outfile = $ARGV[1];
my @files=();

open C, ">$outfile";

my @line = ();

open G, "<$infile";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split "\t", $_;
  	push @files, $line[0];
}

my $a=0;

for ($a=0; $a<@files; $a++) {
	my $sample=((split '/', $files[$a])[1]);
	open F, "<$files[$a]";
	while (<F>){
		chomp;
	  	last if m/^$/;
	  	@line = split "\t", $_;
	  	if ($a == 0) {
	  		if ($line[0] =~ m'AnnotSV') {
		  		print C "Phenotips_ID\t$_\n";
		  	}
	  	}
	  	else {
	  		next if ($line[0] =~ m'AnnotSV');
	  		print C "$sample\t$_\n";
	  	}
	}
	close G;
}