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

# i = integer, s = string
GetOptions ("verbose"           => \$verbose,
            "help"              => \$help
            
      )
 or pod2usage(0) && exit;

pod2usage(1) && exit if ($help);


##############################
# Global variables
##############################


##############################
# Code
##############################


##############################
### Stuff
### More stuff





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

=cut
