#!/usr/bin/perl -w

#takes fasta and writes PRIMER3 input file

use strict;
use Bio::SeqIO;

my $file = shift;
my $inseq = Bio::SeqIO->new(-file   => $file,
                            -format => 'FASTA' );

while (my $seqIO = $inseq->next_seq) {
	my $id  = $seqIO->id;
	my $seq = $seqIO->seq;
	my $len = length($seq);
	#my $optimal_len = $len < 1800 ? $len : 1800;
	my $optimal_len = $len ;
	my $len_75 = int ($optimal_len * .75);
	my $len_50 = int ($optimal_len * .50);
	my $len_25 = int ($optimal_len * .25); 
	print 
"PRIMER_SEQUENCE_ID=$id
SEQUENCE=$seq
PRIMER_PRODUCT_SIZE_RANGE=$len_75-$optimal_len $len_50-$optimal_len $len_25-$optimal_len 100-$len
PRIMER_GC_CLAMP=1	
PRIMER_MIN_GC=40
PRIMER_OPT_GC_PERCENT=50
PRIMER_MAX_GC=60
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60.0
PRIMER_MIN_TM=55.0
PRIMER_MAX_TM=65.0
=\n";
}
