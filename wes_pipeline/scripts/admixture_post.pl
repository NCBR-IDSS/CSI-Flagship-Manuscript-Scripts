#!/usr/bin/perl -w
use strict;
use List::Util 'shuffle';

#INPUT

my $key = $ARGV[0];

@chunks = ('1:1-30000001','1:30000001-60000001','1:60000001-90000001','1:90000001-120000001','1:120000001-150000001','1:150000001-180000001','1:180000001-210000001','1:210000001-240000001','1:240000001-249250621','2:1-30000001','2:30000001-60000001','2:60000001-90000001','2:90000001-120000001','2:120000001-150000001','2:150000001-180000001','2:180000001-210000001','2:210000001-240000001','2:240000001-243199373','3:1-30000001','3:30000001-60000001','3:60000001-90000001','3:90000001-120000001','3:120000001-150000001','3:150000001-180000001','3:180000001-198022430','4:1-30000001','4:30000001-60000001','4:60000001-90000001','4:90000001-120000001','4:120000001-150000001','4:150000001-180000001','4:180000001-191154276','5:1-30000001','5:30000001-60000001','5:60000001-90000001','5:90000001-120000001','5:120000001-150000001','5:150000001-180000001','5:180000001-180915260','6:1-30000001','6:30000001-60000001','6:60000001-90000001','6:90000001-120000001','6:120000001-150000001','6:150000001-171115067','7:1-30000001','7:30000001-60000001','7:60000001-90000001','7:90000001-120000001','7:120000001-150000001','7:150000001-159138663','8:1-30000001','8:30000001-60000001','8:60000001-90000001','8:90000001-120000001','8:120000001-146364022','9:1-30000001','9:30000001-60000001','9:60000001-90000001','9:90000001-120000001','9:120000001-141213431','10:1-30000001','10:30000001-60000001','10:60000001-90000001','10:90000001-120000001','10:120000001-135534747','11:1-30000001','11:30000001-60000001','11:60000001-90000001','11:90000001-120000001','11:120000001-135006516','12:1-30000001','12:30000001-60000001','12:60000001-90000001','12:90000001-120000001','12:120000001-133851895','13:1-30000001','13:30000001-60000001','13:60000001-90000001','13:90000001-115169878','14:1-30000001','14:30000001-60000001','14:60000001-90000001','14:90000001-107349540','15:1-30000001','15:30000001-60000001','15:60000001-90000001','15:90000001-102531392','16:1-30000001','16:30000001-60000001','16:60000001-90000001','16:90000001-90354753','17:1-30000001','17:30000001-60000001','17:60000001-81195210','18:1-30000001','18:30000001-60000001','18:60000001-78077248','19:1-30000001','19:30000001-59128983','20:1-30000001','20:30000001-60000001','20:60000001-63025520','21:1-30000001','21:30000001-48129895','22:1-30000001','22:30000001-51304566','X:1-30000001','X:30000001-60000001','X:60000001-90000001','X:90000001-120000001','X:120000001-150000001','X:150000001-155270560','Y:1-30000001','Y:30000001-59373566','MT:1-16569');
my @line = ();

open H, "<$key";
while (<H>){
	chomp;
  	last if m/^$/;
  	@line = split;
	if ($line[0] ne '#') {
		push @sample, $line[1];
	}
}

my $a=0;
my $step=$ARGV[1];
my $batch=1;
my $b=0;
my $counter=1;

for ($a=0; $a<@chunks; $a++) {
	for ($b=0; $b<@sample; $b++) {
		if ($counter==1) {
			my $list = batch . $batch . '_' . $chunks[$a] . '.list';
			open C, ">$list"
		}
		if ($b<$step) {
			my $file = 'gVCFs/' . $sample[$b] . '.' . $chunks[$a] . '.g.vcf.gz';
			print C "$file\n"
		}
		else {
			$counter=1;
		}
}


open G, "<$ancestry";
while (<G>){
	chomp;
  	last if m/^$/;
  	@line = split;
	if ($line[0] ne '#') {
		my $a=0;
		my $status=0;
		for ($a=0; $a<@sample; $a++) {
			if ($name[$step] eq $sample[$a]) {
				$status++;
				last;
			}
		}
		if ($status == 0) {
			my $e=0;
			print C "$name[$step]";
			for ($e=0; $e<@line; $e++) {
				print C "\t$line[$e]";
			}
			print C "\n";
		}
	}
	$step++;
}