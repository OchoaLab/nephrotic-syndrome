# This script finds the intersection of loci on two files by comparing chr:pos data only.
# It loads the smaller of two BIM files into memory (just chr:pos), then scans the second (larger) file and outputs the lines in that file that are match chr:pos loci in the first file.

# NOTE: this analysis is intended to be for merging a microarray dataset with a WGS dataset like 1000 Genomes.
# the assumptions are different for merging two WGS datasets, this script isn't the best idea for that case!

use strict;

# process inputs
my ($file_in_1, $file_in_2, $file_out) = @ARGV;
die "Usage: perl -w $0 <file1.bim> <file2.bim> <output.bim>\n" unless $file_in_1 && $file_in_2 && $file_out;
# make sure extensions are as expected
die "Error: File must have `.bim` extension: $file_in_1\n" unless $file_in_1 =~ /\.bim$/;
die "Error: File must have `.bim` extension: $file_in_2\n" unless $file_in_2 =~ /\.bim$/;
die "Error: File must have `.bim` extension: $file_out\n" unless $file_out =~ /\.bim$/;

my %data1;

# start reading the first file, load key data onto memory
open(my $file_handle_1, '<', $file_in_1) || die "Could not open for reading $file_in_1: $!";
print "Reading file 1: $file_in_1\n";
while ( <$file_handle_1> ) {
    chomp;
    my @data = split /\t/;
    # make sure there are only 6 fields per line
    die "Error: Line `$_` does not have exactly 6 fields.  File: $file_in_1\n" unless @data == 6;
    
    # get two fields of interest
    my $chr = $data[0];
    my $pos = $data[3];
    
    # toss some data we know we won't be able to map
    next if $chr == 0;
    # sex chromosomes can be mapped in theory, but our 1000 Genomes data doesn't have it, so let's save memory instead
    next if $chr > 22;
    next if $pos == 0;

    # now construct query string, save it in searchable set (hash)
    my $key = $chr . ':' . $pos;
    $data1{ $key } = 1;
}
close $file_handle_1;

# report on progress so far
print "" . ( keys %data1 ) ." loci passed filters\n";

# now start reading the second file, and generating the output file!
open(my $file_handle_2, '<', $file_in_2) || die "Could not open for reading $file_in_2: $!";
open(my $file_handle_out, '>', $file_out) || die "Could not open for writing $file_out: $!";
print "Reading file 2: $file_in_2\n";
print "and writing: $file_out\n";
while ( <$file_handle_2> ) {
    chomp;
    my @data = split /\t/;
    # make sure there are only 6 fields per line
    die "Error: Line `$_` does not have exactly 6 fields.  File: $file_in_2\n" unless @data == 6;
    
    # get two fields of interest, form query key
    my $chr = $data[0];
    my $pos = $data[3];
    my $key = $chr . ':' . $pos;
    
    # if this data appeared in the first file, copy original line (from 2nd file) to output
    print $file_handle_out "$_\n" if $data1{ $key };
}
close $file_handle_2;
close $file_handle_out;

