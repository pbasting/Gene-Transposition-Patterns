#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <utility>

using namespace std;

struct geneInfo{
	string oldLocusTag;
	string locusTag;
	string proteinID;
};

struct syntenyResult{
	string sub_prot;
	string query_prot;
	int moved;
	int moved_adjacent;
	int moved_conserved;
	int conserved_both;
};

struct sortedResults{
	vector<pair<string, string> > not_moved;
	vector<pair<string, string> > moved;
	vector<pair<string, string> > moved_adjacent;
	vector<pair<string, string> > moved_conserved;
	vector<pair<string, string> > conserved_both;
};

struct kegInfo{
	string category;
	vector<string> genes;
};

void parseGenBank(ifstream& genBankFile, vector<geneInfo>& parsedInfo);
void parseKegFile(ifstream& kegFile, vector<kegInfo>& parsedKeg);
void parseSyntenyResults(ifstream& syntenyResults, vector<syntenyResult>& parsedInfo);
void removeMismatches(vector<syntenyResult>& forward, vector<syntenyResult>& reverse);
void findMutualConserved(vector<syntenyResult>& forward, vector<syntenyResult>& reverse);
void sortResults(vector<syntenyResult> forward, vector<syntenyResult> reverse, sortedResults& results);

string parseValue(string line);

int main(int argc, char *argv[]){
	if (argc != 5){
		cout << "missing/too many arguments!"<< endl;
		return 0;
	}
	
	ifstream genBankFile, kegFile ,forwardSyntenyFile, reverseSyntenyFile;
	genBankFile.open(argv[1]);
	kegFile.open(argv[2]);
	forwardSyntenyFile.open(argv[3]);
	reverseSyntenyFile.open(argv[4]);
	
	vector<geneInfo> genBankParsed;
	vector<kegInfo> kegParsed;
	vector<syntenyResult> forwardResults, reverseResults;
	sortedResults resultsSorted;
	
	parseGenBank(genBankFile, genBankParsed);
	parseKegFile(kegFile, kegParsed);
	
	parseSyntenyResults(forwardSyntenyFile, forwardResults);
	parseSyntenyResults(reverseSyntenyFile, reverseResults);
	removeMismatches(forwardResults, reverseResults);
	findMutualConserved(forwardResults, reverseResults);
	sortResults(forwardResults, reverseResults, resultsSorted);

	
	cout << "UNMOVED: " << resultsSorted.not_moved.size() <<endl;
	cout << "MOVED: " << resultsSorted.moved.size() <<endl;
	cout << "MOVED ALONE: " << resultsSorted.moved_adjacent.size() <<endl;
	cout << "MOVED INTO CONSERVED REGION: " << resultsSorted.moved_conserved.size() <<endl;
	cout << "MOVED INTO/FROM CONSERVED REGION: " << resultsSorted.conserved_both.size() <<endl;
	
	return 0;
}


void parseGenBank(ifstream& genBankFile, vector<geneInfo>& parsedInfo){
	string line;
	geneInfo gene;
	while(!genBankFile.eof()){
		//cout <<line <<endl;
		getline(genBankFile, line);
		if (line.find("    CDS    ")!=string::npos){ //finds CDS
			//cout << line <<endl;
			gene.locusTag = "";
			gene.oldLocusTag = "";
			gene.proteinID="";
			do{
				getline(genBankFile, line);
				if(line.find("/locus_tag=") != string::npos){ //finds locus tag
					gene.locusTag = parseValue(line);
					//cout << line <<endl;
				}
				if(line.find("/old_locus_tag=") != string::npos){ //finds locus tag
					gene.oldLocusTag = parseValue(line);
					//cout << line <<endl;
				}
				if(line.find("/protein_id=") != string::npos){ //find protein ID
					gene.proteinID = parseValue(line);
					gene.proteinID = gene.proteinID.substr(0, gene.proteinID.rfind(".")); //removes version number
					//cout << line <<endl;
				}
			}while((line.find("    gene    ")==string::npos) && (!genBankFile.eof()));
			parsedInfo.push_back(gene); //adds struct to vector
		}
	}
}

string parseValue(string line){
	int pos = line.find("=");
	string value;
	pos+=2; //moves past '=' and '"'
	for (pos; ((line[pos] != '\n') && (line[pos] != '"')); pos++){
		value+=line[pos];
	}
	return value;
}

/////////////////////////////////////////////////////////////////
////////////////FINISH WRITING THIS FUNCTION/////////////////////
/////////////////////////////////////////////////////////////////
void parseKegFile(ifstream& kegFile, vector<kegInfo>& parsedKeg){
	string line;
	bool pastOverview = false;
	while(!kegFile.eof()){
		getline(kegFile, line);
		if (line.find("<b>Overview</b>")==string::npos && line.find("<b>") != string::npos && line[0] == 'B'){
			cout << line <<endl;
			pastOverview = true;
		}
		if (pastOverview && line[0] == 'D'){
			cout << line.substr(0, 13) <<endl;
		}
	}	
}

void parseSyntenyResults(ifstream& syntenyResults, vector<syntenyResult>& parsedInfo){
	string line, temp;
	syntenyResult result;
	int pos;
	while(!syntenyResults.eof()){
		getline(syntenyResults, line);
		temp = "";
		if(line.find("lcl|")!=string::npos){ 
		
			//gets subject protein ID////
			pos = line.find("_prot_");
			pos+=6;
			for (pos; ((line[pos] != '\n') && (line[pos] != ',')); pos++){
				temp += line[pos];
			}
			temp = temp.substr(0,temp.rfind("."));
			//cout << temp <<"\t";
			result.sub_prot = temp;
			
			//gets query protein ID///
			temp = "";
			pos = line.rfind("_prot_");
			pos+=6;
			for (pos; ((line[pos] != '\n') && (line[pos] != ',')); pos++){
				temp += line[pos];
			}
			temp = temp.substr(0,temp.rfind("."));
			//cout << temp <<"\t";
			result.query_prot = temp;
			
			//gets movement info//
			line = line.substr(line.find(temp));
			line = line.substr(line.find(",")+1); //moves past query protein
			line = line.substr(line.find(",")+1); // moves past subject protein pos
			line = line.substr(line.find(",")+1); // moves past query protein pos
			result.moved = line[0]-48; //-48 to convert char to int
			//cout <<line[0] << "\t";
			line = line.substr(line.find(",")+1); //moves past 'moved'
			result.moved_adjacent = line[0]-48; //-48 to convert char to int
			//cout <<line[0] << "\t";
			line = line.substr(line.find(",")+1); //moves past 'moved adjacent'
			result.moved_conserved = line[0]-48; //-48 to convert char to int
			//cout <<line[0] << endl;
			result.conserved_both = 0;
			parsedInfo.push_back(result);
		}
	}
}

void removeMismatches(vector<syntenyResult>& forward, vector<syntenyResult>& reverse){
	int count = 0;
	for (int i =0; i < forward.size(); i++){
		for (int j = 0; j < reverse.size(); j++){
			if (forward[i].sub_prot == reverse[j].query_prot){
				if (forward[i].query_prot != reverse[j].sub_prot){
					forward.erase(forward.begin()+i);
					reverse.erase(reverse.begin()+j);
					count++;
				}else if((forward[i].moved != reverse[j].moved) || (forward[i].moved_adjacent != reverse[j].moved_adjacent)){
					forward.erase(forward.begin()+i);
					reverse.erase(reverse.begin()+j);
					count++;
				}
			}
		}
	}
	cout << "mismatches removed: " << count <<endl;
	cout << "remaining forward: " << forward.size() << endl;
	cout << "remaining reverse: " << reverse.size() << endl;
}

void findMutualConserved(vector<syntenyResult>& forward, vector<syntenyResult>& reverse){
	int count = 0;
	for (int i =0; i < forward.size(); i++){
		for (int j = 0; j < reverse.size(); j++){
			if (forward[i].sub_prot == reverse[j].query_prot){
				if (forward[i].moved ==1 && forward[i].moved_adjacent ==1 &&forward[i].moved_conserved ==1 &&
				    reverse[j].moved ==1 && reverse[j].moved_adjacent ==1 &&reverse[j].moved_conserved ==1){
					
					forward[i].conserved_both = 1;
					reverse[j].conserved_both = 1;
					cout <<forward[i].sub_prot << " " << reverse[j].sub_prot <<endl;
					count++;
				}
			}
		}
	}
	cout << "mutually moved into/from conserved: " << count <<endl;
}

void sortResults(vector<syntenyResult> forward, vector<syntenyResult> reverse, sortedResults& results){
	pair<string, string> temp;
	for (int i = 0; i < forward.size(); i++){
		if (forward[i].moved == 0){
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.not_moved.push_back(temp);
		}
		if (forward[i].moved ==1 && forward[i].moved_adjacent == 0){
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.moved.push_back(temp);
		}
		if (forward[i].moved ==1 && forward[i].moved_adjacent == 1 && forward[i].moved_conserved == 0){
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.moved_adjacent.push_back(temp);
		}
		if (forward[i].moved ==1 && forward[i].moved_adjacent == 1 && forward[i].moved_conserved == 1 && forward[i].conserved_both == 0){
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.moved_conserved.push_back(temp);
		}
		if (forward[i].moved ==1 && forward[i].moved_adjacent == 1 && forward[i].moved_conserved == 1 && forward[i].conserved_both == 1){
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.conserved_both.push_back(temp);
		}
		
		//gets the pairs that moved into a conserved region from the reverse perspective
		if (reverse[i].moved ==1 && reverse[i].moved_adjacent == 1 && reverse[i].moved_conserved == 1 && reverse[i].conserved_both == 0){
			temp.first = reverse[i].query_prot;
			temp.second = reverse[i].sub_prot;
			results.moved_conserved.push_back(temp);
		}
	}
}

