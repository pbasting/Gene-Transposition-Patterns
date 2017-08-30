/***************************************************************************************************
FormatKegResults
Written by: Preston Basting
Lab: Jan Mrazek
Last Changed: 8/30/2017
Purpose: This is a component of a series of programs designed to classify protein
		 'movement' when comparing two organisms and determine if proteins belonging
		 to different functional categories are more likely to 'move'
		 
		 This function takes the concatenated 'getKegResults' results of all of the comparisions
		 in a genus and counts the number of hits for each kegg protein category based upon the movement 
		 category generates a table where each row is a different kegg category and each column is a 
		 different movement category. The values are the number or proteins that meat both classifications.
		 
Arguments: (1)any keg file, (2)concatenated genus results from 'getKegResults' (3) name of output .csv
****************************************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <utility>
#include <algorithm>


using namespace std;

struct data{ //stores count data for keg categories
	vector<string> categories;
	vector<int> notMoved, movedAbsolute, movedAdjacent, movedConserved, mutualConserved;
};

vector<string> parseKegFile(ifstream& kegFile);
string parseKegLine(string line);
string upperCase(string line);
void buildTable(vector<string> categories, ifstream& results, data& countData);
void getCategoryCounts(vector<int>& counts, string movement, vector<string> categories, string line, ifstream& file);
void displayTable(data countData);
void outputTable(data countData, ofstream& outputFile);

int main(int argc, char *argv[]){
	if (argc != 4){
		cout << "missing/too many arguments!"<< endl;
		return 0;
	}
	ifstream kegFile, kegResults;
	kegFile.open(argv[1]);
	kegResults.open(argv[2]);
	
	vector<string> kegCategories = parseKegFile(kegFile);
	kegCategories.push_back("UNCATEGORIZED");
	data countData;
	
	
	buildTable(kegCategories, kegResults, countData);
	displayTable(countData);
	
	ofstream output;
	output.open(argv[3]);
	
	outputTable(countData, output);
	

	
	return 0;
}


vector<string> parseKegFile(ifstream& kegFile){
	string line;
	vector<string> categories;
	while(!kegFile.eof()){
		getline(kegFile, line);
		if(line.find("<b>Overview</b>")==string::npos && line[0] == 'B' && line.find("<b>") != string::npos){
			line = parseKegLine(line);
			categories.push_back(line);
		}
	}
	return categories;	
}

string parseKegLine(string line){
	int pos;
	string parsedLine = "";
	pos = line.find("<b>");
	pos+=3;
	while(line[pos] != '<' && pos < line.length()){
		if(line[pos] != ','){ //prevents ',' in category from messing up CSV
			parsedLine += line[pos];
			pos++;
		}else{
			pos++;
		}
	}
	
	return upperCase(parsedLine);
}

string upperCase(string line){
	transform(line.begin(), line.end(), line.begin(), ::toupper);
	return line;
}

void buildTable(vector<string> categories, ifstream& results, data& countData){
	string line="";
	string assignment;
	vector<int> notMoved, movedAbsolute, movedAdjacent, movedConserved, mutualConserved;
	for(int i = 0; i < categories.size(); i++){
			notMoved.push_back(0);
			movedAbsolute.push_back(0);
			movedAdjacent.push_back(0);
			movedConserved.push_back(0);
			mutualConserved.push_back(0);
	}
	
	while(!results.eof()){
		getline(results, line);
		if(line.find("!!")!=string::npos){
			getCategoryCounts(notMoved, "NOT_MOVED", categories, line, results);
			getCategoryCounts(movedAbsolute, "MOVED_ABSOLUTE", categories, line, results);
			getCategoryCounts(movedAdjacent, "MOVED_ADJACENT", categories, line, results);
			getCategoryCounts(movedConserved, "MOVED_CONSERVED", categories, line, results);
			getCategoryCounts(mutualConserved, "MOVED_MUTUAL_CONSERVED", categories, line, results);
		}else{
			getline(results, line);
		}
	}
	
	countData.categories = categories;
	countData.notMoved = notMoved;
	countData.movedAbsolute = movedAbsolute;
	countData.movedAdjacent = movedAdjacent;
	countData.movedConserved = movedConserved;
	countData.mutualConserved = mutualConserved;
	
}


void getCategoryCounts(vector<int>& counts, string movement, vector<string> categories, string line, ifstream& file){
	if(line.find(movement)!=string::npos){
		getline(file, line);
		while(line.find("!!")==string::npos && !file.eof()){
			getline(file, line);
			for (int x = 0; x < categories.size(); x++){
				if (line.find(categories[x]) != string::npos){
					counts[x]+=1;
				}
			}
		}
	}
}

void displayTable(data countData){
	for (int x = 0; x < countData.categories.size(); x++){
		cout << countData.categories[x] << "\t";
		cout << countData.notMoved[x] << "\t";
		cout << countData.movedAbsolute[x] << "\t";
		cout << countData.movedAdjacent[x] << "\t";
		cout << countData.movedConserved[x] << "\t";
		cout << countData.mutualConserved[x] << endl;
	}
}

void outputTable(data countData, ofstream& outputFile){
	outputFile << "FUNCTION,UNMOVED,MOVED.ABS,MOVED.ADJ,MOVED.CONS,MUTUAL.CONS"<<endl;
	for (int x = 0; x < countData.categories.size(); x++){
		outputFile << countData.categories[x] << ",";
		outputFile << countData.notMoved[x] << ",";
		outputFile << countData.movedAbsolute[x] << ",";
		outputFile << countData.movedAdjacent[x] << ",";
		outputFile << countData.movedConserved[x] << ",";
		outputFile << countData.mutualConserved[x] << endl;
	}
}





