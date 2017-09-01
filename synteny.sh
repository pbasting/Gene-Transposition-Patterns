#!/bin/sh

############################################################################################
#synteny.sh
#Written by: Preston Basting
#Email:pjb68507@uga.edu
#Lab: Jan Mrazek
#Last Changed: 9/1/2017
#Purpose: This is a script made to run a series of programs designed to classify protein
#		 'movement' when comparing two organisms and determine if proteins belonging
#		 to different functional categories are more likely to 'move'
#		 
#		 This script runs all the programs necessary for a single pairwise comparison
#		 between two fasta files. It makes blast databases for both fasta files,
#		 runs blastp in both directions, runs 'CompareOrthologs' to get the synteny results,
#		 Makes a synteny map from these results, and runs 'getKegResults' with the synteny results
#		 to classify the proteins by function
#	
#
#Arguments: (1)Subject .fasta (protein), (2)Query .fasta (protein), (3)subject genbank
#			(4) subject .keg
###########################################################################################


sequence1=$1
sequence2=$2
genbank=$3
keg=$4

db_path="databases/"
results_path="blast_results/"
synteny_path="synteny_results/"


seq1Name=${sequence1%.*} #removes file extension
seq1Name=${seq1Name##*/} #removes path to file
seq2Name=${sequence2%.*}
seq2Name=${seq2Name##*/}

if [ ! -d "$db_path" ]; then
	mkdir $db_path
fi

if [ ! -d "$results_path" ]; then
	mkdir $results_path
fi

if [ ! -d "$synteny_path" ]; then
	mkdir $synteny_path
fi

#makes the blast databases
db_1=${db_path}$seq1Name"/"$seq1Name #sequence 1 database
if [ ! -d "${db_path}$seq1Name" ]; then
	mkdir ${db_path}$seq1Name
	makeblastdb -in $sequence1 -out ${db_1}
else
	echo $seq1Name" database already exists"
fi

db_2=${db_path}$seq2Name"/"$seq2Name #sequence 2 database
if [ ! -d "${db_path}$seq2Name" ]; then
	mkdir ${db_path}$seq2Name
	makeblastdb -in $sequence2 -out ${db_2}
else
	echo $seq2Name" database already exists"
fi



blast_results_dir=${results_path}${seq1Name}"_and_"${seq2Name} #location for storage of blast results
blast_results_1=${blast_results_dir}"/subject_"${seq2Name}"_query_"${seq1Name}".txt"
blast_results_2=${blast_results_dir}"/subject_"${seq1Name}"_query_"${seq2Name}".txt"


if [ ! -d "$blast_results_dir" ]; then
	echo "Running blastp....."
	mkdir $blast_results_dir
	#runs blast both directions
	./runBlastp.sh $sequence1 $db_2 $blast_results_1 &
	P1=$!
	./runBlastp.sh $sequence2 $db_1 $blast_results_2 &
	P2=$!
	wait $P1 $P2
else
	echo "blast already run"
fi
	

synteny_dir=${synteny_path}${seq1Name}"_and_"${seq2Name}
mkdir ${synteny_dir}

echo "Finding Moved Proteins..."
#compares ortholog positions and finds proteins that moved
./CompareOrthologs\
 $sequence1\
 $sequence2\
 $blast_results_1\
 $blast_results_2\
 ${synteny_dir}"/subject_"${seq2Name}"_query_"${seq1Name}"_MovementResults.csv"\
 ${synteny_dir}"/subject_"${seq1Name}"_query_"${seq2Name}"_MovementResults.csv"
 
#makes the synteny charts both directions
echo "Making Synteny Plots..."
Rscript makeSyntenyPlot.r ${synteny_dir}"/" "subject_"${seq2Name}"_query_"${seq1Name}"_MovementResults.csv" "subject_"${seq2Name}"_query_"${seq1Name}"_syntenyMap.pdf"
Rscript makeSyntenyPlot.r ${synteny_dir}"/" "subject_"${seq1Name}"_query_"${seq2Name}"_MovementResults.csv" "subject_"${seq1Name}"_query_"${seq2Name}"_syntenyMap.pdf"

#call program that parses genbank and .keg, and converts CompareOrtholog results to keg categories

echo "Assigning KEGG classifications..."
echo " "
./getKegResults\
 $genbank\
 $keg\
 ${synteny_dir}"/subject_"${seq1Name}"_query_"${seq2Name}"_MovementResults.csv"\
 ${synteny_dir}"/subject_"${seq2Name}"_query_"${seq1Name}"_MovementResults.csv"\
 ${synteny_dir}"/kegCounts.csv"
	


