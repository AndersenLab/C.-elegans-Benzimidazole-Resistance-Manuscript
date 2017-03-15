#!/usr/local/bin/perl 

## THIS IS A SCRIPT TO CONVERT BETWEEN WORMBASE GENE IDENTIFIERS FOR C.ELEGANS GENES
## Damien O'Halloran 2015 George Washington University

use strict;
use warnings;
use Getopt::Std;

my %opts;
getopt( 'ab', \%opts );
my $in_name  = $opts{a};
my $out_name = $opts{b};

my $allnames   = "names.csv";
my $input      = "list.txt";
my $outputfile = "converted.txt";
my %unique;
my $all;


open my $in, '<', $input or die "Can't open $input, Perl says $!\n";
open my $out, '>>', $outputfile or die "Can't open $outputfile, Perl says $!\n";

while ( my $line = <$in> ) {
    my $regex = $line;
    chomp $regex;
    $regex =~ s/\s//g;

    open $all, '<', $allnames or die "Can't open $allnames, Perl says $!\n";
    while ( my $line = <$all> ) {

        if (   $line =~ m/,$regex,(.*),/i
            && $in_name eq "WBID"
            && $out_name eq "gene" )
        {
            print $out "$1\n" unless ( $unique{$1}++ );
            $unique{$1} += 1;
        }
        elsif ($line =~ m/,$regex,.*,(.*),/i
            && $in_name eq "WBID"
            && $out_name eq "transcript" )
        {
            print $out "$1\n" unless ( $unique{$1}++ );
            $unique{$1} += 1;
        }
        elsif (   $line =~ m/,$regex,(.*),/i
            && $in_name eq "gene"
            && $out_name eq "transcript" )
        {
            print $out "$1\n" unless ( $unique{$1}++ );
            $unique{$1} += 1;
        }
        elsif ($line =~ m/.*,(.*),$regex,/i
            && $in_name eq "transcript"
            && $out_name eq "gene" )
        {
            print $out "$1\n" unless ( $unique{$1}++ );
            $unique{$1} += 1;
        }
        elsif ($line =~ m/,(.*),$regex,/i
            && $in_name eq "gene"
            && $out_name eq "WBID" )
        {
            print $out "$1\n" unless ( $unique{$1}++ );
            $unique{$1} += 1;
        }
        elsif ($line =~ m/,(.*),.*,$regex,/i
            && $in_name eq "transcript"
            && $out_name eq "WBID" )
        {
            print $out "$1\n" unless ( $unique{$1}++ );
            $unique{$1} += 1;
        }
    }
}

close $in;
close $out;
close $all;
