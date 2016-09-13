# nGlycosylChecker


## Overview

This program takes a protein amino acid sequence and a text file containing information on previously identified missense mutations as input and compares them to output information on any sites where N-glycosylation is affected by a mutation. The output includes the original mutation information, if the glycosylation site is gained or lost, and the amino acid sequence within and flanking the site. 

Any amino acids not matching between the input files will generate an error message with details. 

## Usage

perl nGlycosylChecker.pl [options] --protein protein.fa --mutations mutations.txt > output.txt


## Options

--verbose
	
- Print out additional information to the screen.
	
--help
	
- Print full usage information.
	
-p/--protein

- Fasta file of single protein amino acid sequence.

-m/--mutations

- File containing the mutations to be tested. The input will be included in the output and the third column has to contain the protein change in the format: "p.Ala111Thr". Three letter amino acid codes must be used. An example input file is available in the exampleInput folder. 


## Output

The output adds two columns to the mutation input file. 
- The first indicated if the N-glycosylation site is gained or lost. 
- The next column provides the N-glycosylation site, along with three amino acids flanking both up- and down-stream. The normal and mutated site are both provided, with the normal sequence on the left hand side of the "->" arrow. 