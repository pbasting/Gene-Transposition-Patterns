#!/bin/sh

g++ CompareOrthologs.cpp -o CompareOrthologs

sequence1=$1
sequence2=$2

db_path="databases/"
results_path="blast_results/"
synteny_path="synteny_results/"


seq1Name=${sequence1%.*} #removes file extension
seq1Name=${seq1Name##*/} #removes path to file
seq2Name=${sequence2%.*}
seq2Name=${seq2Name##*/}

#makes the blast databases
db_1=${db_path}$seq1Name"/"$seq1Name #sequence 1 database
if [ ! -d "${db_path}$seq1Name" ]; then
	mkdir ${db_path}$seq1Name
	makeblastdb -in $sequence1 -out ${db_1}	
fi

db_2=${db_path}$seq2Name"/"$seq2Name #sequence 2 database
if [ ! -d "${db_path}$seq2Name" ]; then
	mkdir ${db_path}$seq2Name
	makeblastdb -in $sequence2 -out ${db_2}
fi



blast_results_dir=${results_path}${seq1Name}"_and_"${seq2Name} #location for storage of blast results
blast_results_1=${blast_results_dir}"/subject_"${seq2Name}"_query_"${seq1Name}".txt"
blast_results_2=${blast_results_dir}"/subject_"${seq1Name}"_query_"${seq2Name}".txt"

###############################################
#run these at the same time to speed things up#
###############################################
if [ ! -d "$blast_results_dir" ]; then
	mkdir $blast_results_dir
	#runs blast both directions
	./runBlastp.sh $sequence1 $db_2 $blast_results_1 &
	P1=$!
	./runBlastp.sh $sequence2 $db_1 $blast_results_2 &
	P2=$!
	wait $P1 $P2
fi
	

synteny_dir=${synteny_path}${seq1Name}"_and_"${seq2Name}
mkdir ${synteny_dir}

#compares ortholog positions and finds proteins that moved
./CompareOrthologs\
 $sequence1\
 $sequence2\
 $blast_results_1\
 $blast_results_2\
 ${synteny_dir}"/subject_"${seq2Name}"_query_"${seq1Name}"_MovementResults.csv"\
 ${synteny_dir}"/subject_"${seq1Name}"_query_"${seq2Name}"_MovementResults.csv"\
 ${synteny_dir}"/MutualMovedProteins.csv" 
 
#makes the synteny charts both directions
Rscript makeSyntenyPlot.r ${synteny_dir}"/" "subject_"${seq2Name}"_query_"${seq1Name}"_MovementResults.csv" "subject_"${seq2Name}"_query_"${seq1Name}"_syntenyMap.pdf"
Rscript makeSyntenyPlot.r ${synteny_dir}"/" "subject_"${seq1Name}"_query_"${seq2Name}"_MovementResults.csv" "subject_"${seq1Name}"_query_"${seq2Name}"_syntenyMap.pdf"

#opens the first chart
#evince ${synteny_dir}"/subject_"${seq2Name}"_query_"${seq1Name}"_syntenyMap.pdf"
	


