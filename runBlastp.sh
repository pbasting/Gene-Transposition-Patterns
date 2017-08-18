#!/bin/sh
sequence=$1
database=$2
out_file=$3

blastp\
 -query $sequence\
 -db $database\
 -evalue 1E-20\
 -out $out_file\
 -outfmt "6 qseqid sseqid evalue pident"\
 -max_target_seqs 1\
 -num_threads 8

