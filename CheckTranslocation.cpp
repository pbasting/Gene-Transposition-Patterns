#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
using namespace std;

struct proteinAlignment{
	string qProtein;
	string sProtein;
	string eValue;
	double percentIdentity; 
} proteinData;

void parseBlastResults(ifstream& results, vector<proteinAlignment>& parsedResults);

string getProteinName(string line);

void displayAllResults(vector<proteinAlignment>& parsedResults);

vector<string> getProteinOrder(ifstream& fasta);

void checkForTranslocation(int proteinsToCheck, double threshold, vector<proteinAlignment> proteinData, vector<string> queryProteins, vector<string> subjectProteins);

//int countMatches(vector<string> nearbyProteinsBlast, vector<string>nearbyProteinsFasta);

const int CHECK_RANGE = 5; //number of upstream/downstream proteins to check for translocation

const double THRESHOLD = 0.6; // cutoff for what is considered to have not moved

////////////////////////////////////MAIN/////////////////////////////////////////////////

int main(int argc, char *argv[]){
	if (argc != 4){
		cout << "missing files! Provide: blast results, query fasta, and subject fasta"<< endl;
		return 0;
	}
	vector<proteinAlignment> parsedResults; //vector of all protein alignment results
	ifstream blastResults, queryFasta, subjectFasta;
	blastResults.open(argv[1]); //opens blast output file
	parseBlastResults(blastResults, parsedResults); //parses data into vector of structs
	//displayAllResults(parsedResults); //displays all parsed results to console
	blastResults.close();
	
	queryFasta.open(argv[2]);
	subjectFasta.open(argv[3]);
	vector<string> queryProteinOrder = getProteinOrder(queryFasta);
	vector<string> subjectProteinOrder = getProteinOrder(subjectFasta);
	
	checkForTranslocation(CHECK_RANGE, THRESHOLD, parsedResults, queryProteinOrder, subjectProteinOrder);
	
	return 0;
}






////////////////////////////////////Functions///////////////////////////////////////////

//parses tab-delimited data and stores in approprate variables
void parseBlastResults(ifstream& results, vector<proteinAlignment>& parsedResults){
	for (string line; getline(results, line);){
		int pos;
		for (int x = 0; x < 4; x++){
			pos = line.find("\t");
			switch(x){	//can add more case statements as needed
				case 0: proteinData.qProtein = getProteinName(line.substr(0, pos));
				case 1: proteinData.sProtein = getProteinName(line.substr(0, pos));
				case 2: proteinData.eValue = line.substr(0, pos);
				case 3: proteinData.percentIdentity = atof(line.substr(0, line.length()).c_str());//converting to double
			}
			line = line.substr(pos+1, line.length()-pos);
		}
		parsedResults.push_back(proteinData);
	}
}

string getProteinName(string line){ //parses out protein IDs from seqIDS
	string proteinName = "";
	int namePos = line.find("_prot_");
	for (int x = (namePos+6); line[x-2] != '.'; x++){
		proteinName+=line[x];
	}
	return proteinName;
}

//outputs all of the data contained in the vector of structs
void displayAllResults(vector<proteinAlignment>& parsedResults){
	for (int x = 0; x < parsedResults.size(); x++){
		cout << parsedResults[x].qProtein << "\t";
		cout << parsedResults[x].sProtein << "\t";
		cout << parsedResults[x].eValue << "\t";
		cout << parsedResults[x].percentIdentity << endl;
	}
}

//uses a .fasta to make a vector of protein IDs in order of position in genome
vector<string> getProteinOrder(ifstream& fasta){
	vector<string> proteinOrder;
	for (string line; getline(fasta, line);){
		if (line.find(">") != string::npos){ //finds seqID
			proteinOrder.push_back(getProteinName(line)); //adds protein name to vector of proteins in order of location
			//cout << getProteinName(line)<< endl;
		}
	}
	return proteinOrder;
}

void checkForTranslocation(int proteinsToCheck, double threshold, vector<proteinAlignment> proteinData, vector<string> queryProteins, vector<string> subjectProteins){
	vector<string> qBlastNearbyProteins, qFastaNearbyProteins; //nearby proteins in the blast results and the fasta for the query protein
	//for (int x = 0; x < proteinData.size(); x++){ //uncomment for actual running
	int x = 0;
		for (int y = 0; y < queryProteins.size(); y++){
			if (proteinData[x].qProtein == queryProteins[y]){
				//cout << queryProteins[y] << endl;
				int matches = 0;
				for (int z =1; z <= proteinsToCheck; z++){
					int qDataPos = x;
					int qfastaPos = y;
					qBlastNearbyProteins.push_back(proteinData[x+z].qProtein);
					qFastaNearbyProteins.push_back(queryProteins[y+z]);
					if (x-z >= 0){
						qBlastNearbyProteins.insert(qBlastNearbyProteins.begin(), proteinData[x-z].qProtein); //puts at front of vector
					} else{
						qBlastNearbyProteins.insert(qBlastNearbyProteins.begin(), proteinData[(proteinData.size()) + (x-z)].qProtein); //negative direction in circular chromosome
					}
					if (y-z >= 0){
						qFastaNearbyProteins.insert(qFastaNearbyProteins.begin(), queryProteins[y-z]); //puts at front of vector
					} else{
						qFastaNearbyProteins.insert(qFastaNearbyProteins.begin(), queryProteins[(queryProteins.size()) + (y-z)]);
					}
					
				}
			}
		}
		//the next few lines can be erased
		cout << "BLAST" << "\t\t" << "FASTA" << endl;
		for (int w = 0; w < 10; w++){
			cout << qBlastNearbyProteins[w] << "\t";
			cout << qFastaNearbyProteins[w] << endl;
		}
	//}
}






