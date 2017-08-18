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
	string qPosition;
	string sPosition;
	bool hasMoved; 
} proteinData;

void parseBlastResults(ifstream& results, vector<proteinAlignment>& parsedResults);

string getProteinName(string line);

void parseProteinLocations(string line, string proteinName, vector<proteinAlignment>& parsedResults); //adds subject and query positions to struct

void displayAllResults(vector<proteinAlignment>& parsedResults);

vector<string> getProteinOrder(ifstream& fasta, vector<proteinAlignment>& parsedResults);

void checkForTranslocation(int proteinsToCheck, double threshold, vector<proteinAlignment> proteinData, vector<string> queryProteins, vector<string> subjectProteins, vector<proteinAlignment>& parsedResults);

bool isTranslocated(double threshold, vector<string> nearbyBlastProteins, vector<string> nearbyFastaProteins);

void outputToFile(vector<proteinAlignment>& parsedResults);

void getStats(vector<proteinAlignment>& parsedResults, vector<string> queryProteins, vector<string> subjectProteins);

const int CHECK_RANGE = 6; //number of upstream/downstream proteins to check for translocation

const double THRESHOLD = 0.5; // cutoff for what is considered to have not moved



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
	
	blastResults.close();
	
	queryFasta.open(argv[2]);
	subjectFasta.open(argv[3]);
	vector<string> queryProteinOrder = getProteinOrder(queryFasta, parsedResults); //parses out a list of proteins in order of position from a fasta
	vector<string> subjectProteinOrder = getProteinOrder(subjectFasta, parsedResults);
	
	checkForTranslocation(CHECK_RANGE, THRESHOLD, parsedResults, queryProteinOrder, subjectProteinOrder, parsedResults); //checks all subject proteins for evidence of movement in genome
	displayAllResults(parsedResults); //displays all parsed results to console
	outputToFile(parsedResults);
	getStats(parsedResults, queryProteinOrder, subjectProteinOrder);
	
	
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
		cout << parsedResults[x].percentIdentity << "\t";
		cout << parsedResults[x].hasMoved << "\t";
		cout << parsedResults[x].qPosition << "\t";
		cout << parsedResults[x].sPosition << endl;
	}
}

//uses a .fasta to make a vector of protein IDs in order of position in genome
vector<string> getProteinOrder(ifstream& fasta, vector<proteinAlignment>& parsedResults){
	vector<string> proteinOrder;
	for (string line; getline(fasta, line);){
		if (line.find(">") != string::npos){ //finds seqID
			proteinOrder.push_back(getProteinName(line)); //adds protein name to vector of proteins in order of location
			//cout << getProteinName(line)<< endl;
			parseProteinLocations(line, getProteinName(line), parsedResults);
			
		}
	}
	return proteinOrder;
}

void parseProteinLocations(string line, string proteinName, vector<proteinAlignment>& parsedResults){
	int pos = line.find("[location=");
	if (pos != string::npos){
		pos+=10;
		int tempPos = pos;
		for (int x = 0; x < parsedResults.size(); x++){
			pos = tempPos;
			if (parsedResults[x].qProtein == proteinName){
				while(line[pos-1] != ')' && line[pos] != ']'){
					parsedResults[x].qPosition.push_back(line[pos]);
					pos++;
				}
			}
			else if (parsedResults[x].sProtein == proteinName){
				while(line[pos-1] != ')' && line[pos] != ']'){
					parsedResults[x].sPosition.push_back(line[pos]);
					pos++;
				}
			}
			
		}
	}
}


void checkForTranslocation(int proteinsToCheck, double threshold, vector<proteinAlignment> proteinData, vector<string> queryProteins, vector<string> subjectProteins, vector<proteinAlignment>& parsedResults){
	vector<string> sBlastNearbyProteins, sFastaNearbyProteins; //nearby proteins in the blast results and the fasta for the query protein
	for (int x = 0; x < parsedResults.size(); x++){ //uncomment for actual running
		//int x = 0;
		sBlastNearbyProteins.clear();
		sFastaNearbyProteins.clear();
		for (int y = 0; y < subjectProteins.size(); y++){
			if (parsedResults[x].sProtein == subjectProteins[y]){
				//cout << queryProteins[y] << endl;
				int matches = 0;
				for (int z =1; z <= proteinsToCheck; z++){
					int sDataPos = x;
					int sfastaPos = y;
					//the following finds the porteins downstream of the match
					if (x+z >= parsedResults.size()){
						sBlastNearbyProteins.push_back(proteinData[x-(parsedResults.size())+z].sProtein); //connects end of chromosome w/ start
					}else{
						sBlastNearbyProteins.push_back(parsedResults[x+z].sProtein);
					}
					if (y+z >= subjectProteins.size()){
						sFastaNearbyProteins.push_back(subjectProteins[y-(subjectProteins.size())+z]); //connects end of chroosome w/ start
					}else{
						sFastaNearbyProteins.push_back(subjectProteins[y+z]);
					}
					
					//the following finds the proteins upstream of the match
					if (x-z >= 0){
						sBlastNearbyProteins.insert(sBlastNearbyProteins.begin(), parsedResults[x-z].sProtein); //puts at front of vector
					}else{
						sBlastNearbyProteins.insert(sBlastNearbyProteins.begin(), parsedResults[(proteinData.size()) + (x-z)].sProtein); //negative direction in circular chromosome
					}
					if (y-z >= 0){
						sFastaNearbyProteins.insert(sFastaNearbyProteins.begin(), subjectProteins[y-z]); //puts at front of vector
					}else{
						sFastaNearbyProteins.insert(sFastaNearbyProteins.begin(), subjectProteins[(subjectProteins.size()) + (y-z)]);
					}
					
				}
			}
		}
		//the next few lines can be erased
		/*cout << "******" << parsedResults[x].sProtein << "******" << endl;
		cout << "BLAST" << "\t\t" << "FASTA" << endl;
		for (int w = 0; w < 10; w++){
			cout << sBlastNearbyProteins[w] << "\t";
			cout << sFastaNearbyProteins[w] << endl;
		}
		cout << isTranslocated(threshold, sBlastNearbyProteins, sFastaNearbyProteins) << endl;
		*/
		parsedResults[x].hasMoved = isTranslocated(threshold, sBlastNearbyProteins, sFastaNearbyProteins);
	}
}

bool isTranslocated(double threshold, vector<string> nearbyBlastProteins, vector<string> nearbyFastaProteins){
	double matches = 0.0;
	for (int x = 0; x < nearbyBlastProteins.size(); x++){
		for (int y = 0; y < nearbyFastaProteins.size(); y++){
			if (nearbyBlastProteins[x] == nearbyFastaProteins[y]){
				matches++;
			}
		}
	}
	if (matches/double(nearbyBlastProteins.size()) < threshold){
		return true;
	} else {
		return false;
	}
}

void outputToFile(vector<proteinAlignment>& parsedResults){
	ofstream outputFile;
	outputFile.open ("results.txt");
	outputFile << "query ID" << "\t" << "subject ID" << "\t" << "evalue" << "\t" << "Percent Identity" << "\t";
	outputFile << "Moved" << "\t" << "query location" << "\t" << "subject location" << endl;
	for (int x = 0; x < parsedResults.size(); x++){
		outputFile << parsedResults[x].qProtein << "\t";
		outputFile << parsedResults[x].sProtein << "\t";
		outputFile << parsedResults[x].eValue << "\t";
		outputFile << parsedResults[x].percentIdentity << "\t";
		if (parsedResults[x].hasMoved == true){
			outputFile << "true" << "\t";
		}else{
			outputFile << "false" << "\t";
		}
		outputFile << parsedResults[x].qPosition << "\t";
		outputFile << parsedResults[x].sPosition << endl;
	}
}

void getStats(vector<proteinAlignment>& parsedResults, vector<string> queryProteins, vector<string> subjectProteins){
	cout << parsedResults.size() <<"/" << queryProteins.size() << " aligned" << endl;
	int count = 0;
	for (int x = 0; x < parsedResults.size(); x++){
		if (parsedResults[x].hasMoved == true){
			count++;
		}
	}
	cout << count << " were predicted to have moved in the genome" << endl;
}





