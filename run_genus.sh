#!/bin/sh

#the only argument passed is the name of the genus
#this should match the directory where the fastas are stored in fastas/
genus=$1
g=${genus:0:1}

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
	cat ${dir}/kegCounts.csv >> synteny_results/${genus}/${genus}_movedProteins.txt #concatenates all results for a genus into a .csv
	mv ${dir} synteny_results/${genus}
done

#moves all the blast results to a single directory named after the genus
mkdir blast_results/${genus}
mv blast_results/${g}_* blast_results/${genus}


./FormatKegResults $keg1 synteny_results/${genus}/${genus}_movedProteins.txt synteny_results/${genus}/${genus}_formatted_movedProteins.csv
