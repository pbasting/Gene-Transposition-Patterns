#!/bin/sh

#the only argument passed is the name of the genus
#this should match the directory where the fastas are stored in fastas/
genus=$1
g=${genus:0:1}

#runs the synteny script with every pair of fastas for a genus
for file1 in fastas/${genus}/${g}_*.fasta; do
	for file2 in fastas/${genus}/${g}_*.fasta; do
	file1Name=${file1%.*} #removes file extension
	file1Name=${file1Name##*/} #removes path to file
	file2Name=${file2%.*}
	file2Name=${file2Name##*/}
		if [ $file1 != $file2 ] && [ ! -d "blast_results/"$file2Name"_and_"$file1Name ] #makes sure this hasen't been run before
		then
			echo $file1Name and $file2Name
			./synteny.sh $file1 $file2
		fi 
	done
done

mkdir synteny_results/${genus}
touch synteny_results/${genus}/${genus}_movedProteins.csv

#moves all synteny_results into a directory named after the genus
for dir in synteny_results/${g}_*; do
	cat ${dir}/MutualMovedProteins.csv >> synteny_results/${genus}/${genus}_movedProteins.csv #concatenates all results for a genus into a .csv
	mv ${dir} synteny_results/${genus}
done

#moves all the blast results to a single directory named after the genus
mkdir blast_results/${genus}
mv blast_results/${g}_* blast_results/${genus}
