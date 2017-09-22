#!/usr/bin/python

############################################################################################
#ConstructBrKegg.py
#Written by: Preston Basting
#Email:pjb68507@uga.edu
#Lab: Jan Mrazek
#Last Changed: 9/22/2017
#Purpose: This is a script made to construct a reference file of genes sorted into functional
#		  categories. It takes a kegg orthology file and downloads the appropriate BRITE file
#		  It constructs a single file using the information from both
#	
#Arguments: (1) .keg file (2) output file name (.brkeg)
###########################################################################################
import urllib2
import sys



##gets kegg file
inputFile = open(sys.argv[1], "r")
inputLines = inputFile.readlines()
inputFile.close()

#finds keg org name
for x in range (0,len(inputLines)):
	if (inputLines[x].find("#ENTRY") != -1):
		pos = inputLines[x].find("#ENTRY")
		pos+=7
		org = ''
		while pos < len(inputLines[x]):
			if (inputLines[x][pos].isalpha()):
				org+=inputLines[x][pos]
			pos+=1
			

print("Downloading "+org+" KEGG results...")
#gets brite file from the web
brURL = 'http://rest.kegg.jp/link/'+org+'/brite'
brite = urllib2.urlopen(brURL, timeout=20).read()

#parses brite file
line = ''
briteLines = []
while x < len(brite):
	while brite[x] != "\n":
		line += brite[x]
		x+=1
	briteLines.append(line)
	line = ''
	x+=1

#stores brite file in dictionary. key = br category values = genes
from collections import defaultdict
briteDict = defaultdict(list)
for x in range (0, len(briteLines)):
	br = briteLines[x][briteLines[x].find(":")+1:]
	gene = br[br.find(":")+1:]
	br = br[:br.find("\t")]
	briteDict[br].append(gene)


#makes dictionary of all kegg categories and the associated genes including the brite genes
geneList = []
keggDict = defaultdict(list)
for x in range(0, len(inputLines)):
	if (inputLines[x].find("</b>") !=-1 and inputLines[x][0] == 'B' and inputLines[x].find("Overview") == -1):
		category = inputLines[x][:inputLines[x].find("</b>")]
		category = category[category.find("<b>")+3:]
		y = x+1
		while inputLines[y][0] != 'B' and inputLines[y][0] !='!':
			if inputLines[y][0] == 'C' and inputLines[y].find("[BR:") != -1:
				br = inputLines[y][inputLines[y].find(":")+1:]
				br = br[:br.find("]")]
				if br in briteDict:
					for q in briteDict[br]:
						keggDict[category].append(q)
			if inputLines[y][0] == 'D':
				gene = inputLines[y][inputLines[y].find(" ")+6:]
				gene = gene[:gene.find(" ")]
				keggDict[category].append(gene)
			y+=1
		


print("Building .brkeg file...")
#outputs the data to a file for downstream parsing and use
output = open(sys.argv[2], "w")
for categories in keggDict:
	output.write("C\t<")
	output.write(categories)
	output.write(">\n")
	for genes in keggDict[categories]:
		output.write("G\t<")
		output.write(genes)
		#output.write(" ")
		output.write(">\n")
			
output.close()

print("Done with "+org)

