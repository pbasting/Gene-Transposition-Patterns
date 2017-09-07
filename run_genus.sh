#!/bin/sh

############################################################################################
#run_genus.sh
#Written by: Preston Basting
#Email:pjb68507@uga.edu
#Lab: Jan Mrazek
#Last Changed: 9/7/2017
#Purpose: This is a script designed to run another script made to run a series of programs 
#		  designed to classify protein 'movement' when comparing two organisms and determine 
#		  if proteins belonging to different functional categories are more likely to 'move'
#		 
#		 This script is used to run 'synteny.sh' for each subdirectory of a given directory.
#		 It also concatenates the 'getKegResults' output into a single file for the entire
#		 directory. This is then parsed and formatted by 'FormatKegResults' to produce a
#		 .csv file containing a table of keg category counts vs movement category.
#
#Arguments: (1)Name of directory/genus (no path, must be within 'fastas' directory)
###########################################################################################


g++ CompareOrthologs.cpp -o CompareOrthologs
g++ getKegResults.cpp -o getKegResults
g++ FormatKegResults.cpp -o FormatKegResults

#the only argument passed is the name of the genus
#this should match the directory where the fastas are stored in fastas/
genus=$1
g=${genus:0:1}

#runs avery pairwise comparison without duplicates
for org1 in fastas/${genus}/${g}_*; do
	fasta1=$org1/*.fasta
	gb1=$org1/*.gb
	keg1=$org1/*.keg
	for org2 in fastas/${genus}/${g}_*; do
		fasta2=$org2/*.fasta
		gb2=$org2/*.gb
		keg2=$org2/*.keg
		file1Name=${org1%.*} #removes file extension
		file1Name=${file1Name##*/} #removes path to file
		file2Name=${org2%.*}
		file2Name=${file2Name##*/}
		if [ $org1 != $org2 ] && [ ! -d "blast_results/"$file2Name"_and_"$file1Name ] #makes sure this hasen't been run before
		then
			echo $file1Name and $file2Name
			#echo $fasta1 and $fasta2
			./synteny.sh $fasta1 $fasta2 $gb1 $keg1
		fi
	done
done



mkdir synteny_results/${genus}
touch synteny_results/${genus}/${genus}_movedProteins.txt

#moves all synteny_results into a directory named after the genus
for dir in synteny_results/${g}_*; do
	cat ${dir}/kegCounts.csv >> synteny_results/${genus}/${genus}_movedProteins.txt #concatenates all results for a genus into a single file
	mv ${dir} synteny_results/${genus}
done

#moves all the blast results to a single directory named after the genus
mkdir blast_results/${genus}
mv blast_results/${g}_* blast_results/${genus}

#parses the concatenated results, counts the hits for each keg category base upon movement category, and stores in a .csv as a table
echo "Formatting genus KEGG results..."
./FormatKegResults $keg1 synteny_results/${genus}/${genus}_movedProteins.txt synteny_results/${genus}/${genus}_formatted_movedProteins.csv

#uses poisson distribution to determine probability of category counts occuring
echo "Running Poisson approximations..."
Rscript getPoissonValues.r synteny_results/${genus}/${genus}_formatted_movedProteins.csv synteny_results/${genus}/${genus}_poisson.csv synteny_results/${genus}/${genus}_Summary_Table.csv


