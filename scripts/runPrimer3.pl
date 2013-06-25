#!/usr/bin/perl -w

# this script will write a primer3 input file and will run parsePrimer3.pl 
# runPrimer3.pl -f fasta.nt  (Note -f is required)
# input a fasta file
# output a tab delimited file will all possible primers

use strict;
use Getopt::Std;
use Bio::SeqIO;


#####################
#
# reading, storing cmd line variables
#
###################


my %options=();
my $optString = 'h:g:t:l:f:';
sub init(){
        getopts($optString,\%options) or usage();
        # like the shell getopt, "d:" means d takes an argument

        usage() if defined $options{h};
        usage() if !defined $options{f};
	if (defined $options{g} and $options{g} !~ /\D+/) {
                usage();
        }
	if (defined $options{t} and $options{t} !~ /\D+/) {
                usage();
        }
	if (defined $options{l} and $options{l} !~ /\D+/) {
                usage();
        }
	
        print "Unprocessed by Getopt::Std:\n" if $ARGV[0];
}

sub usage() {
        print STDERR <<EOF;

To run this program  ...

usage: $0 [-h] [-g GC-opt,min,max] [-t Tm-opt,min,max][-l PrimerLen-opt,min,max][-f file]

required:
        -f file : fasta file
optional:
        -h      : this (help) message
        -g      : GC content -- optimal,min,max (default is 50,45,60)
        -t      : Melting Temperature -- optimal,min,max (default is 55,50,65)
        -l      : Primer Length  -- optimal min max (default is 20,18,25)

        example: $0 -f seqs.fasta
        example: $0 -g 50,45,60 -t 60,55,65 -f seqs.fasta


EOF
	exit;
}

init();


####################
#
# tests 
# 
####################

my $file = $options{f};

#test to see if file exists
die "\n\n***Oops the file: \'$file\' does not exist.\n Is it spelled right??\n\n" unless -e $file;

#replace carriage returns for newlines
`perl -p -i -e 's/\r/\n/g' $file`;

#test for fasta format
my $found = `grep ">" $file`;

if (!$found){
	die "\n\n***ERROR***
		Input file not in FASTA format.
		No primers were produced.

	EXAMPLE FASTA FILE

	>SEQNAME
	GTATATATGGGTAAATAAATTGAAAATTACTTGTCTTAT
	TCCTTCAGTAATTGCCGTAACTGAAAATTTTGATGCTGA
	ATAGAAATGAATTGAAGAAAAGCTGAAAACCTTGTGACC	


"
;
}
####################
#
#  stores command line variables into program variables
#    -- or --
#  if no cmd line variables are given it stores default values
#
####################


my ($gcOpt, $gcMin, $gcMax) = defined $options{g} ? split /\D+/, $options{g} : (50,45,60) ;
my ($tmOpt, $tmMin, $tmMax)  = defined $options{t} ? split /\D+/, $options{t} : (55,50,65) ;
my ($primerLenOpt, $primerLenMin, $primerLenMax) = defined $options{l} ? split /\D+/, $options{l} : (20,18,25) ;


####################
#
#  1.parsing of the required fasta file
#  2.writing the primer3 input file
#
####################

my $inseq = Bio::SeqIO->new(-file   => $file,
                            -format => 'FASTA' );
my $count = 0;
while (my $seqIO = $inseq->next_seq) {
        my $id  = $seqIO->id;
        my $seq = $seqIO->seq;
        my $len = length($seq);
        my $len_75 = int ($len * .75);
        my $len_50 = int ($len * .50);
        my $len_25 = int ($len * .25);

	
	if ($len_25 <= $primerLenMax){
		warn "\n\n  ***** WARNING ****** $id sequence is too short.  No primers will be made for this sequence.\n\n";
		next if $len_25 <= $primerLenMax;
	}	
	

$file =~ s/\.txt//;
open OUTFILE,  ">>$file.in_primer3";

        print OUTFILE

"PRIMER_SEQUENCE_ID=$id
SEQUENCE=$seq
PRIMER_PRODUCT_SIZE_RANGE=$len_75-$len $len_50-$len $len_25-$len 100-$len
PRIMER_GC_CLAMP=1       
PRIMER_MIN_GC=$gcMin
PRIMER_OPT_GC_PERCENT=$gcOpt
PRIMER_MAX_GC=$gcMax
PRIMER_OPT_SIZE=$primerLenOpt
PRIMER_MIN_SIZE=$primerLenMin
PRIMER_MAX_SIZE=$primerLenMax
PRIMER_OPT_TM=$tmOpt
PRIMER_MIN_TM=$tmMin
PRIMER_MAX_TM=$tmMax
=\n";
	$count++
}



if ($count){ #proceed only if 1 or more appropriate sequences are provided

####################
#
# running the primer3
#
###################

	print "Primer3 is finding primers for the *$count* seqeunces in $file\n\n";

	`cat $file.in_primer3  | primer3_core > $file.out_primer3`;

	print "Your Primers are in $file.primers.txt\n\n";



####################
#
# running the parser
#
###################

	print `/usr/lib/cgi-bin/primer3/parsePrimer3.pl $file.out_primer3 > $file.primers.txt`;
}
