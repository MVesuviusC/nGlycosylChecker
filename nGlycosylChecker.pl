#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;

##############################
# By Matt Cannon
# Date: 08/25/16
# Last modified: 08/25/16
# Title: nGlycosylChecker.pl
# Purpose: 
##############################

##############################
# Options
##############################


my $verbose;
my $help;
my $protein;
my $mutations;

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help,
	    "protein=s"         => \$protein,
	    "mutations=s"       => \$mutations        
	    )
    or pod2usage(0) && exit;

# print out full usage if asked
pod2usage(1) && exit if ($help);
# print out short usage if two required files not provided
pod2usage(0) && exit if(!defined($protein) || !defined($mutations));

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
}
close PROTEIN;

##############################
### Pull in list of mutations and check for changed N-glyc. sites

local $/ = "\n";

open MUTATIONS, $mutations or die "Cannot open mutation list file\n";
while(my $line = <MUTATIONS>) {
    chomp $line;
    my ($ntChange, $genPos, $aaChange) = split "\t", $line;
    $aaChange =~ s/^p.//;
    my $pos = substr($aaChange, 3);
    $pos =~ s/...$//;
    #print STDERR $aaChange, "\t", $pos, "\t";
    my $origAA = lc(substr($aaChange, 0, 3));
    #print STDERR $origAA, "\t";
    my $mutAA = lc(substr($aaChange, -3, 3));
    #print STDERR $mutAA, "\t";
    my $origSeq = substr($sequence, $pos - 3, 5);
    #print STDERR $origSeq, "\t";
    my $mutSeq = substr($sequence, $pos - 3, 2) . $aaCodeHash{$mutAA} . substr($sequence, $pos, 2);
    #print STDERR $mutSeq, "\n";
    if($origSeq =~ /N.[XT]/ && $mutSeq !~ /N.[XT]/) {
	print $line, "\tLost_N-glycosylation_site\n";
    } elsif($origSeq !~ /N.[XT]/ && $mutSeq =~ /N.[XT]/) {
	print $line, "\tGained_N-glycosylation_site\n";
    }

    
}



##############################
# POD
##############################

#=pod
    
=head SYNOPSIS

Summary:    
    
    nGlycosylChecker.pl - Checks a mutation to see if it adds or removes an N-glycosylation site
    
Usage:

    perl nGlycosylChecker.pl [options] 


=head OPTIONS

Options:

    --verbose
    --help
    --protein
    --mutations

=cut
