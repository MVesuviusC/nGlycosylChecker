#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

##############################
# By Matt Cannon
# Date: 08/25/16
# Last modified: 09/13/16
# Title: nGlycosylChecker.pl
# Purpose: Find mutations that change N-glycosylation sites
##############################

##############################
# Options
##############################

my $verbose;
my $help;
my $protein;
my $mutations;

GetOptions (
	    "verbose"           => \$verbose,
            "help"              => \$help,
	    "protein=s"         => \$protein,
	    "mutations=s"       => \$mutations        
	    )
    or pod2usage(0) && exit;

# print out full usage if asked
pod2usage(-verbose => 1) && exit if ($help);
# print out short usage if two required files not provided
my $msg = "Please provide both required files\n\n";
pod2usage(-verbose => 0, -msg => $msg) && exit if(!defined($protein) || !defined($mutations));

##############################
# Global variables
##############################
my $singleSeqCheck = 0;
my $sequence;
my %aaCodeHash = (
		  ala => "A",
		  cys => "C",
		  asp => "D",
		  glu => "E",
		  phe => "F",
		  gly => "G",
		  his => "H",
		  ile => "I",
		  lys => "K",
		  leu => "L",
		  met => "M",
		  asn => "N",
		  pro => "P",
		  gln => "Q",
		  arg => "R",
		  ser => "S",
		  thr => "T",
		  val => "V",
		  trp => "W",
		  tyr => "Y"
		  );
my $nGlyRegex = "N[ACDEFGHIKLMNQRSTVWY][ST]";

##############################
# Code
##############################


##############################
### Pull in protein sequence
# Change line delimiter to pull in entire fasta entry
local $/ = ">\n";
open PROTEIN, $protein or die "Cannot open protein fasta file\n";
while(my $line = <PROTEIN>) {
    chomp $line;
    # Make sure there is only one protein sequence
    if($singleSeqCheck == 1 && $line =~ /^>g/) {
	print STDERR "More than one fasta entry in protein fasta file\n";
	print STDERR "Please provide only one entry\n";
	die;
    }

    my ($header, @rawSeq) = split "\n", $line;
    $sequence = join("", @rawSeq);

    ###########
    ### QC
    # Kick out any spaces or tabs in the sequence
    $sequence =~ s/[\s\t]//g;
    # Make sure they're all uppercase
    $sequence = uc($sequence);

    # If any unknown amino acids are in the sequence, die
    if($sequence =~ /[BJOUXZ]/) {
	print STDERR "Illegal amino acid symbol in protein sequence [BJOUXZ]\n\n";
	die;
    }
    ###
    ###########
    $singleSeqCheck = 1;

    # Add 5 spaces to the beginning and end of $sequence
    # To avoid issues where the mutation is at the beginning
    # or end of the sequence
    $sequence = "     " . $sequence . "     ";


}
close PROTEIN;

##############################
### Pull in list of mutations and check for changed N-glyc. sites

local $/ = "\n";

open MUTATIONS, $mutations or die "Cannot open mutation list file\n";
while(my $line = <MUTATIONS>) {
    chomp $line;
    my (undef, undef, $aaChange) = split "\t", $line;
    $aaChange =~ s/^p.//;
    
    ## Check if $aaChange is premature stop
    if($aaChange =~ /\*/) {
	print $line, "\n\n\nInvalid mutation type (premature stop)!\n\n\n";
	die;
    }

    my $pos = substr($aaChange, 3);
    $pos =~ s/...$//;
    # Make sure position is a number
    if($pos =~ /[a-zA-Z]/) {
	print STDERR "\n\n\nError in mutation information!\n", $line, "\n\n\n";
	die;
    }

    # Add 5 to $pos to make up for five spaces in front of $sequence
    $pos += 5;

    # get original and mutant proteins from mutation info
    my $origAA = lc(substr($aaChange, 0, 3));
    my $mutAA = lc(substr($aaChange, -3, 3));

    # Check that $origAA and $mutAA are amino acids
    if(!defined($aaCodeHash{$origAA}) || !defined($aaCodeHash{$mutAA})) {
	print STDERR "\n\n\nError in mutation information!\n", $line, "\n\n\n";
	die;
    }
    
    # Get five five AA sequence surrounding mutated AA (2 on either side)
    #    - Unless the AA is within 3 of the ends of the sequence, in which case
    #      grab all sequences within 2 of the mutant (but don't grab from the 
    #      other end).
    my $origSeq;
    if($pos <= 3) {
	$origSeq = substr($sequence, 0, $pos + 2);
	#print STDERR $line, "\t", $origSeq, "\t", $sequence, "\n";
    } elsif ((length($sequence) - $pos) <= 3 ) {
	$origSeq = substr($sequence, $pos - 3, length($sequence) - $pos);
	#print STDERR $line, "\t", $origSeq, "\t", $sequence, "\n";
    } else {
	$origSeq = substr($sequence, $pos - 3, 5);
    }

    # Double check that the protein seq. and coordinates line up
    my $checkAA = substr($sequence, $pos - 1, 1);
    if( $checkAA ne $aaCodeHash{$origAA}) {
	print STDERR 
	    $line, "\t", $origSeq, "\t.", $sequence, ".\n",
	    "Amino acid in fasta sequence \"", 
	    $checkAA, 
	    "\" does not match amino acid in mutation file \"", 
	    $aaCodeHash{$origAA}, "\"\n\n";
		
    }


    my $mutSeq = substr($sequence, $pos - 3, 2) . $aaCodeHash{$mutAA} . substr($sequence, $pos, 2);

    my $upStreamSeq = substr($sequence, $pos - 8, 5);
    my $downStreamSeq = substr($sequence, $pos + 2, 5);

    if($verbose) {
	print STDERR $line, "\t", $origSeq, " -> ", $mutSeq, "\n";
    }
    
    # Check if glycosylation site changes and print out any that do
    if($origSeq =~ /$nGlyRegex/ && $mutSeq !~ /$nGlyRegex/) {
	my $origInFront = substr($origSeq, 0, $-[0]);
	my $origMatch = substr($origSeq, $-[0], 3);
	my $mutMatch = substr($mutSeq, $-[0], 3);
	my $origInBack = substr($origSeq, $+[0], length($origSeq) - $+[0]);
	my $origCombined = 
	    substr($upStreamSeq . $origInFront, -3, 3) . 
	    "_" . $origMatch . "_" . 
	    substr($origInBack . $downStreamSeq, 0, 3);
	my $mutCombined = 
	    substr($upStreamSeq . $origInFront, -3, 3) .
	    "_" . $mutMatch . "_" . 
	    substr($origInBack . $downStreamSeq, 0, 3);
	print $line, "\tLost_N-glycosylation_site\t", $origCombined, "->", $mutCombined, "\n";
    } elsif($origSeq !~ /$nGlyRegex/ && $mutSeq =~ /$nGlyRegex/) {
	my $origInFront = substr($origSeq, 0, $-[0]);
	my $origMatch = substr($origSeq, $-[0], 3);
	my $mutMatch = substr($mutSeq, $-[0], 3);
	my $origInBack = substr($origSeq, $+[0], length($origSeq) - $+[0]);
	my $origCombined = 
	    substr($upStreamSeq . $origInFront, -3, 3) . 
	    "_" . $origMatch . "_" . 
	    substr($origInBack . $downStreamSeq, 0, 3);
	my $mutCombined = 
	    substr($upStreamSeq . $origInFront, -3, 3) .
	    "_" . $mutMatch . "_" . 
	    substr($origInBack . $downStreamSeq, 0, 3);
	print $line, "\tGained_N-glycosylation_site\t", $origCombined, "->", $mutCombined, "\n";
    }
}
close MUTATIONS;


##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    nGlycosylChecker.pl - Checks a list of mutations to see if they add or remove N-glycosylation sites
    
Usage:

    perl nGlycosylChecker.pl [options] --protein protein.fa --mutations mutations.txt > output.txt


=head OPTIONS

Options:

    --verbose

        Print out additional information to the screen.

    --help
    
        Print full usage information.

    -p/--protein
    
        Fasta file of single protein sequence.

    -m/--mutations

        File containing the mutations to be tested. The input will be included in the output
        and the third column has to contain the protein change in the format: "p.Ala111Thr". Three letter 
        amino acid codes must be used. 

=cut
