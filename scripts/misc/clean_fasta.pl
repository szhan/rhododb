#! /usr/bin/perl
use strict;
use warnings;

use lib "/data/home/szhan/projects/misc/1kp-algae/scripts/";
use CommonTools qw(parse_fasta);


my $USAGE = "USAGE: perl $0 <in fasta> <out fasta>\n";
die $USAGE unless scalar @ARGV == 2;

my $in_fa_file = $ARGV[0];
my $out_fa_file = $ARGV[1];


my %SEQS = %{ &parse_fasta($in_fa_file) };

open OUT, ">$out_fa_file" or die "ERROR: failed to create $out_fa_file";
foreach ( keys %SEQS ) {
	my $id = $_;
	my $seq = $SEQS{$_};

	$seq =~ s/\-//g;
	print OUT ">$id\n$seq\n";
}
close OUT;

