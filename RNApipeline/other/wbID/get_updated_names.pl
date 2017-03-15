#!/usr/bin/perl

use strict;
use warnings;
use LWP::Simple;

my $out_dir = 'WormBase_IDs';
if ( !-e $out_dir ) {
    print "Output directory ($out_dir) does not exist! Creating...\n";
    mkdir($out_dir) or die "Failed to create the output directory ($out_dir): $!\n";
    print "Done\n\n";
}

my $content;
my $out_file;
my $out;

$out_file = 'names.csv';

#find out which version of wormbase you want to use
print
"\nHello, please enter the version of Wormbase you want IDs from - example would be WS239\n";
my $version = <STDIN>;
# Remove the newline
chomp $version;
$version =~ s/\s//;
$version =~ tr/a-z/A-Z/;

print
"Loading and retrieving the list of Celegans IDs from Wormbase $version...\n";
$content = get(
'ftp://ftp.wormbase.org/pub/wormbase/species/c_elegans/annotation/geneIDs/c_elegans.PRJNA13758.'
      . $version
      . '.geneIDs.txt.gz
'
);
die "Failed to download gene names!" unless $content;
open ($out, '>', "$out_dir/$out_file")  or die "Can't create the gene name file ($out_dir/$out_file): $!\n";

print $out "$content\n";

print "\n\nall done....\n\n";

close $out;


