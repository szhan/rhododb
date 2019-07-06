#! /usr/bin/perl
use strict;
use warnings;

my $USAGE = "USAGE: perl $0 <in fasta> <out phylip>\n";
die $USAGE unless scalar @ARGV == 2;

my $in_fa_file = $ARGV[0];
my $out_phy_file = $ARGV[1];

# Print to outfile in Phylip format
open OUTFILE, ">$out_phy_file" or die "ERROR: failed to create $out_phy_file";
print OUTFILE &convert_fasta2phylip($in_fa_file, $out_phy_file);
close OUTFILE;

# Convert MSA FastA to Phylip format
sub convert_fasta2phylip {
	my $in_file = $_[0];
	my $out_file = $_[1];
	
	my $nbr_seq = 0;
	my $nbr_pos = 0;
	my $phy_txt = '';
	
	my $tmp_id = '';
	my $tmp_seq = '';
	
	open INFILE, $in_file or die "ERROR: failed to open $in_file";
	while(<INFILE>){
		s/[\r\n]//g;
		if ( m/^>/ ) {
			# Update
			$nbr_seq++;
			$nbr_pos = length($tmp_seq);
			$phy_txt .= "$tmp_id $tmp_seq\n" if $tmp_id ne '';
			
			# Reset
			$tmp_id = substr($_, 1);
			$tmp_seq = '';
		} else {
			$tmp_seq .= $_;
		}
	}
	close INFILE;
	
	# Add last entry
	# Same as in the Update section
	$phy_txt .= "$tmp_id\t$tmp_seq\n";
	
	# Prepend number of sequences and alignment positions
	$phy_txt = " $nbr_seq $nbr_pos\n$phy_txt";
	
	return $phy_txt;
}

