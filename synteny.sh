#!/bin/sh
sequence1=$1
sequence2=$2

db_path="/home/preston/Documents/blastp_project/databases/"
results_path="/home/preston/Documents/blastp_project/blast_results/"
synteny_path="/home/preston/Documents/blastp_project/synteny_results/"


seq1Name=${sequence1%.*}
seq1Name=${seq1Name##*/}
seq2Name=${sequence2%.*}
seq2Name=${seq2Name##*/}

mkdir ${db_path}$seq1Name
mkdir ${db_path}$seq2Name

makeblastdb -in $sequence1 -out ${db_path}$seq1Name"/"$seq1Name
makeblastdb -in $sequence2 -out ${db_path}$seq2Name"/"$seq2Name

mkdir ${results_path}${seq1Name}"_and_"${seq2Name}

./runBlastp.sh $sequence1 ${db_path}$seq2Name"/"$seq2Name ${results_path}${seq1Name}"_and_"${seq2Name}"/subject_"${seq2Name}"_query_"${seq1Name}".txt"
./runBlastp.sh $sequence2 ${db_path}$seq1Name"/"$seq1Name ${results_path}${seq1Name}"_and_"${seq2Name}"/subject_"${seq1Name}"_query_"${seq2Name}".txt"


synteny_dir=${synteny_path}${seq1Name}"_and_"${seq2Name}
mkdir ${synteny_dir}
./CompareOrthologs\
 $sequence1\
 $sequence2\
 ${results_path}${seq1Name}"_and_"${seq2Name}"/subject_"${seq2Name}"_query_"${seq1Name}".txt"\
 ${results_path}${seq1Name}"_and_"${seq2Name}"/subject_"${seq1Name}"_query_"${seq2Name}".txt"\
 ${synteny_dir}"/subject_"${seq2Name}"_query_"${seq1Name}"_MovementResults.csv"\
 ${synteny_dir}"/subject_"${seq1Name}"_query_"${seq2Name}"_MovementResults.csv"
 
Rscript makeSyntenyPlot.r ${synteny_dir}"/" "subject_"${seq2Name}"_query_"${seq1Name}"_MovementResults.csv" "subject_"${seq2Name}"_query_"${seq1Name}"_syntenyMap.pdf"

Rscript makeSyntenyPlot.r ${synteny_dir}"/" "subject_"${seq1Name}"_query_"${seq2Name}"_MovementResults.csv" "subject_"${seq1Name}"_query_"${seq2Name}"_syntenyMap.pdf"

evince ${synteny_dir}"/subject_"${seq2Name}"_query_"${seq1Name}"_syntenyMap.pdf"
	


