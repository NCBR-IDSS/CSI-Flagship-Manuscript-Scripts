#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $infile = $ARGV[0];
my $id = $ARGV[1];
my @files=();
my $cmd='';

my @line = ();

open G, "<$infile";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split ',', $_;
  	if ($line[0] ne 'error') {
		if ($line[7] eq $id) {
			if ($line[6] eq 'male') {
				$cmd = "cp CSI_wes_pipeline/resources/male_ploidy.vcf CANVAS/$id/ploidy.vcf";
				system($cmd);
			}
			else {
				$cmd = "cp CSI_wes_pipeline/resources/female_ploidy.vcf CANVAS/$id/ploidy.vcf";
				system($cmd);
			}
		}
	}
}