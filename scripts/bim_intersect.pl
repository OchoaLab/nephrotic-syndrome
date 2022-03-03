# This script finds the intersection of loci on two files by comparing chr:pos data only minimally (can also compare alleles in case where they are expected to be aligned).
# It assumes loci are ordered in both, reads them together and outputs common lines to two different outputs (to account for fact that IDs, alleles, etc, may differ)
# NOTE: assumes no multiallelics (no two rows in a single file with same chr/pos), if this is not true, some matches will be accidentally excluded

use strict;

# a constant
# if false, alleles are ignored (only chr/pos are compared)
my $match_alleles = 1;

# process inputs
my ($file_in_1, $file_in_2) = @ARGV;
die "Usage: perl -w $0 <file1.bim> <file2.bim>\nOutputs get suffix '_out' added before extension\n" unless $file_in_1 && $file_in_2;
# make sure extensions are as expected
die "Error: File must have `.bim` extension: $file_in_1\n" unless $file_in_1 =~ /\.bim$/;
die "Error: File must have `.bim` extension: $file_in_2\n" unless $file_in_2 =~ /\.bim$/;

# construct output paths
my $file_out_1 = $file_in_1; # copy
$file_out_1 =~ s/\.bim$//; # remove extension
$file_out_1 .= '_out.bim'; # concatenate suffix and add extension back
my $file_out_2 = $file_in_2; # copy
$file_out_2 =~ s/\.bim$//; # remove extension
$file_out_2 .= '_out.bim'; # concatenate suffix and add extension back

# open input files
open(my $handle_in_1, '<', $file_in_1) || die "Could not open for reading $file_in_1: $!";
open(my $handle_in_2, '<', $file_in_2) || die "Could not open for reading $file_in_2: $!";
# now open outputs
open(my $handle_out_1, '>', $file_out_1) || die "Could not open for writing $file_out_1: $!";
open(my $handle_out_2, '>', $file_out_2) || die "Could not open for writing $file_out_2: $!";

# parse a single line
sub parse_line {
    my ($handle_in, $file) = @_;
    # read line
    my $line = <$handle_in>;
    # stop if file is done
    return ($line) unless defined $line;
    # process line otherwise
    chomp $line;
    my @data = split /\t/, $line;
    # make sure there are only 6 fields per line
    die "Error: Line `$line` does not have exactly 6 fields.  File: $file\n" unless @data == 6;
    # get fields of interest
    my $chr = $data[0];
    my $pos = $data[3];
    my $ref = $data[4];
    my $alt = $data[5];
    # return only the bits we want
    return ($line, $chr, $pos, $ref, $alt);
}

# parse bits
my ($line1, $chr1, $pos1, $ref1, $alt1) = parse_line( $handle_in_1, $file_in_1 );
my ($line2, $chr2, $pos2, $ref2, $alt2) = parse_line( $handle_in_2, $file_in_2 );

# for optional allele checks
my $good = 1;

while ( $line1 && $line2 ) {
    # check for equality
    if ( $chr1 == $chr2 && $pos1 == $pos2 ) {
	# decide whether to print to output
	if ( $match_alleles ) {
	    # print only if both alleles match
	    $good = ( $ref1 eq $ref2 && $alt1 eq $alt2 );
	} else {
	    $good = 1; # print regardless, not checking alleles
	}
	if ( $good ) {
	    # print these lines to their outputs
	    print $handle_out_1 $line1."\n";
	    print $handle_out_2 $line2."\n";
	}
	# either way (whether we printed or not), advance both lines
	($line1, $chr1, $pos1, $ref1, $alt1) = parse_line( $handle_in_1, $file_in_1 );
	($line2, $chr2, $pos2, $ref2, $alt2) = parse_line( $handle_in_2, $file_in_2 );
    } elsif ( $chr1 < $chr2 || ($chr1 == $chr2 && $pos1 < $pos2 ) ) {
	# file 1 is behind, so advance that one
	($line1, $chr1, $pos1, $ref1, $alt1) = parse_line( $handle_in_1, $file_in_1 );
    } else {
	#} elsif ( $chr1 > $chr2 || ($chr1 == $chr2 && $pos1 > $pos2 ) ) {
	# the only other possibility is that file 2 is behind
	($line2, $chr2, $pos2, $ref2, $alt2) = parse_line( $handle_in_2, $file_in_2 );
    }
}

# if we're here, at least one file ended, so we're done!
close $handle_in_1;
close $handle_in_2;
close $handle_out_1;
close $handle_out_2;
