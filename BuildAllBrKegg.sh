#!/bin/sh
############################################################################################
#BuildAllBrKegg.sh
#Written by: Preston Basting
#Email:pjb68507@uga.edu
#Lab: Jan Mrazek
#Last Changed: 9/22/2017
#Purpose: This is a script made to construct the functional databases used downstream to
#		  categorize genes into functions. It calls 'ConstructBrKegg.py' for all keg files
#
#Arguments: None
###########################################################################################
for species in fastas/*/*; do
		speciesName=${species%.*}
		speciesName=${speciesName##*/}
		./ConstructBrKegg.py ${species}/*.keg ${species}/${speciesName}.brkeg
done
	
