#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

vector<string> getProteinIDs(ifstream& fasta);

vector<int> getMatchPositions(ifstream& blast, vector<string> subjectFastaProteins, vector<string> queryFastaProteins);

vector<string> getBlastSubjectProteins(ifstream& blast);

vector<string> getBlastqueryProteins(ifstream& blast);

int findBestMatch(ifstream& blast, vector<int> matchPositions,vector<string> queryFastaProteins);

double getPercentIdentity(string line);

void outputResults( vector<int> matchPositions, vector<string> subjectFasta, vector<string> queryFasta, ofstream& outputFile);

bool checkForMovement(vector<int> matchPositions, int index);

bool checkAdjacentProteins(vector<int> matchPositions, int index, int maxQuerySize);

const int PROTEIN_POS_CUTOFF = 200; //the difference from original protein position that can be considered a large move
const int CHECK_RANGE = 5; //number of upstream and downstream proteins to check for differences
const int RANGE_CUTOFF = 30; //divergence from the checked protein that can still be considered a nearby protein

const double NEARBY_PROTEIN_CUTOFF = 0.4; //cutoff for fraction of different nearby proteins for a protein that hasn't moved //lower = more conservative
const double PERCENT_IDENTITY_CUTOFF = 30.0; //lowest acceptable percent identity for matches
const int NO_PROTEIN = -1; //indicates no protein match in vectors of match positions


int main(int argc, char *argv[]){
	if (argc != 7){
		cout << "missing/too many arguments! Provide:  query fasta, subject fasta, forward blast results, reverse blast results, and output file name"<< endl;
		return 0;
	}
	
	//builds vectors of proteins for query and subject fasta
	ifstream queryFasta, subjectFasta, forwardBlast, reverseBlast;
	ofstream outputFile1, outputFile2;
	queryFasta.open(argv[1]);
	subjectFasta.open(argv[2]);
	forwardBlast.open(argv[3]);
	reverseBlast.open(argv[4]);
	outputFile1.open(argv[5]);
	outputFile2.open(argv[6]);
	if(!queryFasta.is_open() || !subjectFasta.is_open() || !forwardBlast.is_open() || !reverseBlast.is_open() || !outputFile1.is_open() || !outputFile2.is_open()){
		cout << "ERROR:failed to open one of the files" << endl;
		return 0;
	}
	
	vector<string> queryFastaProteins, subjectFastaProteins;
	queryFastaProteins = getProteinIDs(queryFasta);
	subjectFastaProteins = getProteinIDs(subjectFasta);
	
	//builds vectors of positions: index= position of protein in subject fasta, value= position of protein in query fasta
	vector<int> forwardMatchPositions = getMatchPositions(forwardBlast, subjectFastaProteins, queryFastaProteins);
	vector<int> reverseMatchPositions = getMatchPositions(reverseBlast, queryFastaProteins, subjectFastaProteins);
	
	//checks for movement and outputs results
	outputResults(forwardMatchPositions, subjectFastaProteins, queryFastaProteins, outputFile1);
	outputResults(reverseMatchPositions, queryFastaProteins, subjectFastaProteins, outputFile2);	


	return 0;
}

vector<string> getProteinIDs(ifstream& fasta){
	vector<string> proteinIDs;
	for (string line; getline(fasta, line);){
		if (line[0] == '>'){
			int x =1;
			string temp = "";
			while((line[x] != ' ') && (line[x] != '\t') && (x < line.length())){
				temp+=line[x];
				x++;
			}
			//cout << temp << endl;
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
	for (string line; getline(blast, line);){
		temp = "";
		int pos = line.find("\t");
		pos++; //beginging of subject proteinID
		for (pos; line[pos] != '\t' && pos < line.length(); pos++){
			temp+= line[pos];
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
			if (subjectFastaProteins[x] == blastSubjectProteins[y]){
				for(int z=0; z < queryFastaProteins.size(); z++){
					if(blastQueryProteins[y] == queryFastaProteins[z]){
						matches.push_back(z);
					}
				}
			}
		}
		if (matches.size() ==0){
			matchPositions.push_back(NO_PROTEIN);
		} else if (matches.size() ==1){
			matchPositions.push_back(matches[0]);
		} else {
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
	while (line[pos] < line.length()){
		pIdent+=line[pos];
		pos++;
	}
	return atof(pIdent.c_str());
	
	
}

//writes important information to a file comma-delimited
void outputResults( vector<int> matchPositions, vector<string> subjectFasta, vector<string> queryFasta, ofstream& outputFile){
	outputFile << "Subject.Protein" <<","<<"Query.Protein" <<","<< "Movement.Distance" <<","<<"Movement.Adjacent" << endl;
	for (int x = 0; x < matchPositions.size(); x++){

		if (matchPositions[x] >= 0){
			//outputFile << subjectFasta[x] << "\t";
			outputFile << x << ",";
			//outputFile << queryFasta[matchPositions[x]] << "\t";
			outputFile << matchPositions[x] << ",";
			outputFile << checkForMovement(matchPositions, x) << ",";
			outputFile << checkAdjacentProteins(matchPositions,x, queryFasta.size());
		}
		outputFile << endl;
	}
}

///////////////////////////////////////////////////////////////////////////////////
/////think about how to account for circular chromosome when the query and subject are different sizes
///////////////////////////////////////////////////////////////////////////////////
bool checkForMovement(vector<int> matchPositions, int index){
	if (abs(index - matchPositions[index]) > PROTEIN_POS_CUTOFF){
		return true; //protein moved
	}else{
		return false; //protein did not move
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







