#!/usr/bin/perl
use strict;
use warnings;

# Description
# -----------
# Code to reformat CAST output for the covid-19-Annotations-on-Structures
# project https://github.com/gtauriello/covid-19-Annotations-on-Structures

#
# Author: Vasilis J Promponas 
# Contact: vprobon@ucy.ac.cy; vasilis.promponas@gmail.com
#

my %colorscheme = getColorScheme();

my @annot = ();
open(CASTTAB, $ARGV[0]) or die "Could not open CAST tab ($ARGV[0])";
my $entry = <CASTTAB>; # Read header line

while($entry = <CASTTAB>)
{
# Input format (tab separated):
# >sp|P0DTD3|Y14_WCPV Uncharacterized protein 14 OS=Wuhan seafood market pneumonia virus OX=2697049 GN=ORF14 PE=3 SV=1    C       68      70      27
	chomp($entry);
	next if $entry eq '';
	my @tmp = split(/\t/,$entry);
	my $aatype = $tmp[1]. '-rich region (' . $tmp[4]. ')';
	my $from = $tmp[2];
	my $to  = $tmp[3];
	my @tmp2 = split(/\|/, $tmp[0]);
	my $uniprot_id = $tmp2[1];
	my $color = $colorscheme{$tmp[1]};
	my $out = join("\t", ($uniprot_id,  $from, $to, $color, $aatype,)) . "\n";
	push @annot, $out;
}
close(CASTTAB);


foreach my $annotation (@annot)
{
	print $annotation;
}


sub getColorScheme
{
	# Return a hash mapping the Rasmol color scheme
	# See http://life.nthu.edu.tw/~fmhsu/rasframe/SHAPELY.HTM
	my %colors=();
	$colors{'D'} = $colors{'E'} = '#E60A0A';
	$colors{'C'} = $colors{'M'} = '#E6E600';
	$colors{'K'} = $colors{'R'} = '#145AFF';
	$colors{'S'} = $colors{'T'} = '#FA9600';
	$colors{'F'} = $colors{'Y'} = '#3232AA';
	$colors{'N'} = $colors{'Q'} = '#00DCDC';
	$colors{'G'} = '#EBEBEB';
	$colors{'L'} = $colors{'V'} = $colors{'I'} = '#0F820F';
	$colors{'A'} = '#C8C8C8';
	$colors{'W'} = '#B45AB4';
	$colors{'H'} = '#8282D2';
	$colors{'P'} = '#DC9682';
	
	return(%colors);
}
