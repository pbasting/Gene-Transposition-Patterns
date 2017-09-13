/***************************************************************************************************
CompareOrthologs
Written by: Preston Basting
Email:pjb68507@uga.edu
Lab: Jan Mrazek
Last Changed: 9/1/2017
Purpose: This is a component of a series of programs designed to classify protein
		 'movement' when comparing two organisms and determine if proteins belonging
		 to different functional categories are more likely to 'move'
		 
		 This function takes format 6 results from blastp (containing qseqid sseqid evalue pident)
		 done in both directions for a set of organisms. It then classifies the relative movement
		 of the proteins(didn't move, moved with adjacent proteins, moved without adjacent proteins, 
		 moved without adjacent proteins into a conserved region). these are output into two csv files
		 (one for each direction of the blastp) which are used by downstream programs.
		 
Arguments: (1)query.fasta, (2)subject.fasta, (3)forwardBlast, (4)reverseBlast, (5)forwardOutputFile,
		   (6)reverseOutputFile
****************************************************************************************************/
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <utility>

using namespace std;

//Takes a fasta file and generates a vector of all of the proteinIDs
vector<string> getProteinIDs(ifstream& fasta);

//generates a vector of positions where the index corresponds to the protein in the  subject fasta and the value 
//corresponds to the protein in the query fasta. These are linked based upon the blast results
vector<int> getMatchPositions(ifstream& blast, vector<string> subjectFastaProteins, vector<string> queryFastaProteins);

//parses out all of the subject proteinIDs from the blast results and stores in a vector
vector<string> getBlastSubjectProteins(ifstream& blast);

//parses out all of the query proteinIDs from the blast results and stores in a vector
vector<string> getBlastqueryProteins(ifstream& blast);

//if multiple matches exist for the query proteins, the position of the match with the best percent Identity is returned
int findBestMatch(ifstream& blast, vector<int> matchPositions,vector<string> queryFastaProteins);

//parses out the percent Identity of a given match from the blast results
double getPercentIdentity(string line);

//outputs all the protein matches and the movement classification information to a csv
void outputAllResults( vector<int> matchPositions, vector<string> subjectFasta, vector<string> queryFasta, ofstream& outputFile);

//checks upstream and downstream query proteins to see if they are conserved compared to the subject
//returns true if the surrounding proteins are different
bool checkAdjacentProteins(vector<int> matchPositions, int index, int maxQuerySize);

//checks upstream and downstream query proteins to see if the protein entered a conserved region
//returns true if the region is conserved
bool isConserved(vector<int> matchPositions, int index, int maxQuerySize);

//parses out the file name without the path
string getFileName(string fileAndPath);


const int CHECK_RANGE = 5; //number of upstream and downstream proteins to check for differences
const int RANGE_CUTOFF = 5; //divergence from the checked protein that can still be considered a nearby protein

const double NEARBY_PROTEIN_CUTOFF = 0.3; //cutoff for fraction of different nearby proteins for a protein that hasn't moved //lower = more conservative
const double PERCENT_IDENTITY_CUTOFF = 50.0; //lowest acceptable percent identity for matches
const int NO_PROTEIN = -1; //indicates no protein match in vectors of match positions



////////////////////////////MAIN//////////////////////////////////////////////////////
int main(int argc, char *argv[]){
	if (argc != 7){
		cout << "missing/too many arguments! Provide:  query fasta, subject fasta, forward blast results, reverse blast results, and output file name"<< endl;
		return 0;
	}
	
	//builds vectors of proteins for query and subject fasta
	ifstream queryFasta, subjectFasta, forwardBlast, reverseBlast;
	ofstream outputFile1, outputFile2, mutualMatchesOut;
	queryFasta.open(argv[1]);
	subjectFasta.open(argv[2]);
	forwardBlast.open(argv[3]);
	reverseBlast.open(argv[4]);
	outputFile1.open(argv[5]);
	outputFile2.open(argv[6]);
	if(!queryFasta.is_open() || !subjectFasta.is_open() || !forwardBlast.is_open() || !reverseBlast.is_open() || !outputFile1.is_open() || !outputFile2.is_open()){
		cout << "C++ERROR:failed to open one of the files" << endl;
		return 0;
	}
	
	vector<string> queryFastaProteins, subjectFastaProteins;
	queryFastaProteins = getProteinIDs(queryFasta);
	subjectFastaProteins = getProteinIDs(subjectFasta);
	
	//builds vectors of positions: index= position of protein in subject fasta, value= position of protein in query fasta
	vector<int> forwardMatchPositions = getMatchPositions(forwardBlast, subjectFastaProteins, queryFastaProteins);
	vector<int> reverseMatchPositions = getMatchPositions(reverseBlast, queryFastaProteins, subjectFastaProteins);

	//checks for movement and outputs results
	outputAllResults(forwardMatchPositions, subjectFastaProteins, queryFastaProteins, outputFile1);
	outputAllResults(reverseMatchPositions, queryFastaProteins, subjectFastaProteins, outputFile2);
	
	//outputs the proteins that were found to move based on absolute position and adjacent proteins
	//only outputs the proteins that were found in both directions of the blast
	string subjectName = getFileName(argv[2]);
	string queryName = getFileName(argv[1]);	

	return 0;
}


////////////////////FUNCTIONS////////////////////////////////////////////////////////////////

vector<string> getProteinIDs(ifstream& fasta){
	vector<string> proteinIDs;
	for (string line; getline(fasta, line);){
		if (line[0] == '>'){ //finds the line where the annotations are
			int x =1;
			string temp = "";
			while((line[x] != ' ') && (line[x] != '\t') && (x < line.length())){
				temp+=line[x];
				x++;
			}
			proteinIDs.push_back(temp);
		}
	}
	return proteinIDs;
}

vector<string> getBlastSubjectProteins(ifstream& blast){
	vector<string> subjectProteins;
	string temp;
	blast.clear();
	blast.seekg(0, ios::beg);
	int pos;
	for (string line; getline(blast, line);){
		temp = "";
		pos = line.find("\t");
		pos++; //beginging of subject proteinID
		while((line[pos] != '\t') && (line[pos] != ' ') && (pos < line.length())){
			temp+= line[pos];
			pos++;
		}
		if (getPercentIdentity(line) >= PERCENT_IDENTITY_CUTOFF){ //only adds proteins that meet the cutoff
			subjectProteins.push_back(temp);
		}
	}
	return subjectProteins;
}

vector<string> getBlastQueryProteins(ifstream& blast){
	vector<string> queryProteins;
	string temp;
	blast.clear();
	blast.seekg(0, ios::beg);
	for (string line; getline(blast, line);){
		temp="";
		for (int pos = 0; line[pos] != '\t' && pos < line.length(); pos++){
			temp+= line[pos];
		}
		if (getPercentIdentity(line) >= PERCENT_IDENTITY_CUTOFF){ //only adds proteins that meet the cutoff
			queryProteins.push_back(temp);
		}
	}
	return queryProteins;
}

//this function is building a  vector of matches to coorelate subject proteins (indexes) to query proteins (values) based upon fasta positions
vector<int> getMatchPositions(ifstream& blast, vector<string> subjectFastaProteins, vector<string> queryFastaProteins){
	vector<int> matchPositions, matches;
	vector<string> blastQueryProteins, blastSubjectProteins;
	blastSubjectProteins = getBlastSubjectProteins(blast);
	blastQueryProteins = getBlastQueryProteins(blast);
	for (int x = 0; x < subjectFastaProteins.size(); x++){
		matches.clear();
		for (int y =0; y < blastSubjectProteins.size(); y++){
			if (subjectFastaProteins[x] == blastSubjectProteins[y]){ //finds the subject protein in the blast
				for(int z=0; z < queryFastaProteins.size(); z++){
					if(blastQueryProteins[y] == queryFastaProteins[z]){ //finds the query protein in the query fasta
						matches.push_back(z); //adds the position of the query protein at the index of the sujbect protein match
					}
				}
			}
		}
		if (matches.size() ==0){
			matchPositions.push_back(NO_PROTEIN);
		} else if (matches.size() ==1){
			matchPositions.push_back(matches[0]);
		} else { //handles multiple matches//returns the match with the best percent identity
			matchPositions.push_back(findBestMatch(blast, matches, queryFastaProteins));
		}
	}
	return matchPositions;	
}

int findBestMatch(ifstream& blast, vector<int> matchPositions, vector<string> queryFastaProteins){
	int bestMatch;
	double tempPerIdent;
	double topPerIdent = 0; //top percent identity
	for (int x = 0; x < matchPositions.size(); x++){
		blast.clear();
		blast.seekg(0, ios::beg);
		for (string line; getline(blast, line);){
			if (line.find(queryFastaProteins[matchPositions[x]]) != string::npos){
				tempPerIdent = getPercentIdentity(line);
				if (topPerIdent < tempPerIdent){
					topPerIdent = tempPerIdent;
					bestMatch = matchPositions[x];
				}
			}
		}
	}
	return bestMatch;
	
}

double getPercentIdentity(string line){
	int pos = line.rfind("\t");
	string pIdent = "";
	while ((line[pos] < line.length()) && (pIdent.length() <6)){
		pIdent+=line[pos];
		pos++;
	}
	return atof(pIdent.c_str());
	
	
}

//writes important information to a file comma-delimited
void outputAllResults( vector<int> matchPositions, vector<string> subjectFasta, vector<string> queryFasta, ofstream& outputFile){
	outputFile << "S_Prot_Name, Q_Prot_Name,Subject.Protein,Query.Protein,Movement.Adjacent,Adjacent.Conserved"<<endl;
	for (int x = 0; x < matchPositions.size(); x++){

		if (matchPositions[x] >= 0){
			outputFile << subjectFasta[x] << ",";
			outputFile << queryFasta[matchPositions[x]] << ",";
			outputFile << x << ",";
			
			outputFile << matchPositions[x] << ",";
			outputFile << checkAdjacentProteins(matchPositions,x, queryFasta.size()) << ",";
			outputFile << isConserved(matchPositions,x, queryFasta.size());
			outputFile << endl;
		}
		
	}
}


bool checkAdjacentProteins(vector<int> matchPositions, int index, int maxQuerySize){
	int minVal, maxVal;
	double count = 0; //counts proteins that stayed within range
	double totalChecked = (CHECK_RANGE*2);
	vector<int> matchesToCheck;
	
	minVal = matchPositions[index] - RANGE_CUTOFF;
	if (minVal < 1){
		minVal = ((maxQuerySize+matchPositions[index]) - RANGE_CUTOFF);
	}
	
	maxVal = matchPositions[index] + RANGE_CUTOFF;
	if (maxVal > maxQuerySize){
		maxVal = ((matchPositions[index] + RANGE_CUTOFF) - (maxQuerySize));
	}
	
	//makes sub-vector of matches to check//
	
	//gets downstream values
	int i = index+1;
	while (matchesToCheck.size() < CHECK_RANGE){
		if (i < matchPositions.size()){
			if (matchPositions[i] > -1){
				matchesToCheck.push_back(matchPositions[i]);
			}
			i++;
		}else{
			i = 0;
		}
	}
	
	//gets upstream values
	i = index-1;
	while (matchesToCheck.size() < totalChecked){
		if (i >= 0){
			if (matchPositions[i] > -1){
				matchesToCheck.push_back(matchPositions[i]);
			}
			i--;
		}else{
			i = (matchPositions.size()-1);
		}
	}
	
	//checks values in subvector to see if they are adjacet to the same proteins in the subject sequence
	//count is the count of nearby proteins that are the same in both genomes nearby
	if (minVal < maxVal){
		for (int x =0; x < matchesToCheck.size(); x++){
			if (matchesToCheck[x] >= minVal && matchesToCheck[x] <= maxVal){
				count+=1;
			}
		}
	}else{
		for (int x =0; x < matchesToCheck.size(); x++){
			if (matchesToCheck[x] >= minVal || matchesToCheck[x] <= maxVal){
				count+=1;
			}
		}
	}

	
	if ((count/totalChecked) < NEARBY_PROTEIN_CUTOFF){
		return true; //moved
	}else{
		return false; //did not move
	}
}

bool isConserved(vector<int> matchPositions, int index, int maxQuerySize){
	double totalChecked = CHECK_RANGE*2;
	double count = 0;
	//gets upstream proteins//
	int x=index-1;
	vector<int> adjacentProteins;
	while (adjacentProteins.size() < CHECK_RANGE){
		if (x < 0){
			x = matchPositions.size()-1;
		}
		if (matchPositions[x] > -1){
			adjacentProteins.insert(adjacentProteins.begin(), matchPositions[x]);
		}
		x--;
	}
	//gets downstream proteins//
	x = index+1;
	
	while(adjacentProteins.size() < totalChecked){
		if (x >= matchPositions.size()){
			x = 0;
		}
		if(matchPositions[x] > -1){
			adjacentProteins.push_back(matchPositions[x]);
		}
		x++;
	}
	
	vector<int> proteinsToCheck;
	int i =0;
	while (i < maxQuerySize){
		count = 0;
		
		//compares all adjacent proteins to a series of increasing numbers
		//counts the total matches
		//this is looking to see if the adjacent proteins are a conserved region
		int j =i;
		int q = 0;
		while(q < totalChecked*2){
			if (j >= maxQuerySize){
				j = 0;
			}
			for(int k=0; k < adjacentProteins.size(); k++){
				if(j == adjacentProteins[k]){
					count+=1;
				}
			}
			j++;
			q++;
		}
		if((count/totalChecked) > 1.0-NEARBY_PROTEIN_CUTOFF){
			return true;
			
		}
		if (count == 0){
			i+=10; //speeds up the search by skipping areas with no matches
		}else{
			i++;
		}
	}
	return false;
}

string getFileName(string fileAndPath){
	string fileName = "";
	int pos = fileAndPath.rfind("/"); //moves past the file path
	pos++;
	for (pos; fileAndPath[pos] != '.'; pos++){ //doesn't include file extension
		fileName+= fileAndPath[pos];
	}
	return fileName;
}





