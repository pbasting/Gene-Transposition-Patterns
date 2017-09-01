/***************************************************************************************************
FormatKegResults
Written by: Preston Basting
Email:pjb68507@uga.edu
Lab: Jan Mrazek
Last Changed: 9/1/2017
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

//takes a keg file and parses out the categories. Stores all the categories into a vector
vector<string> parseKegFile(ifstream& kegFile);

//takes a line from a keg file taht contains a category. Parses out the name of the category
string parseKegLine(string line);

//converts a string to uppercase
string upperCase(string line);

//takes the results from 'getKegResults', counts how many proteins match to each category and splits
//the count into movement categories
void buildTable(vector<string> categories, ifstream& results, data& countData);

//used by build table to tally counts for each movement category in each keg category
void getCategoryCounts(vector<int>& counts, string movement, vector<string> categories, string line, ifstream& file);

//outputs the count results to a .csv
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
	kegCategories.push_back("UNCATEGORIZED"); //adds uncategorized as a category 
	
	data countData;
	buildTable(kegCategories, kegResults, countData);
	
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
		//skips overview
		//all category lines start with B and are wrapped in <b></b>
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
	pos+=3; //skips past <b>
	while(line[pos] != '<' && pos < line.length()){
		parsedLine += line[pos];
		pos++;
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
	//adds zeros to all vector positions so specific indexes can be increased during the count
	for(int i = 0; i < categories.size(); i++){
			notMoved.push_back(0);
			movedAbsolute.push_back(0);
			movedAdjacent.push_back(0);
			movedConserved.push_back(0);
			mutualConserved.push_back(0);
	}
	
	while(!results.eof()){
		getline(results, line);
		if(line.find("!!")!=string::npos){ //!! indicates the line contains a movement category
			getCategoryCounts(notMoved, "NOT_MOVED", categories, line, results);
			getCategoryCounts(movedAbsolute, "MOVED_ABSOLUTE", categories, line, results);
			getCategoryCounts(movedAdjacent, "MOVED_ADJACENT", categories, line, results);
			getCategoryCounts(movedConserved, "MOVED_CONSERVED", categories, line, results);
			getCategoryCounts(mutualConserved, "MOVED_MUTUAL_CONSERVED", categories, line, results);
		}
	}
	
	//builds countData struct
	countData.categories = categories;
	countData.notMoved = notMoved;
	countData.movedAbsolute = movedAbsolute;
	countData.movedAdjacent = movedAdjacent;
	countData.movedConserved = movedConserved;
	countData.mutualConserved = mutualConserved;
	
}


void getCategoryCounts(vector<int>& counts, string movement, vector<string> categories, string line, ifstream& file){
	if(line.find(movement)!=string::npos){
		while(line.find("**")==string::npos){ //** indicates end of a movement category
			getline(file, line);
			for (int x = 0; x < categories.size(); x++){
				if (line.find(categories[x]) != string::npos){ //if the category is found in the line
					counts[x]+=1;
				}
			}
		}
	}
}


void outputTable(data countData, ofstream& outputFile){
	outputFile << "FUNCTION,UNMOVED,MOVED.ABS,MOVED.ADJ,MOVED.CONS,MUTUAL.CONS"<<endl;
	for (int x = 0; x < countData.categories.size(); x++){
		if(countData.categories[x] == upperCase("Folding, sorting and degradation")){ //the comma in this category messes up the comma delimiting
			outputFile << upperCase("Folding sorting and degradation,"); //this is the same catagory title without the comma
		}else{
			outputFile << countData.categories[x] << ",";
		}
		outputFile << countData.notMoved[x] << ",";
		outputFile << countData.movedAbsolute[x] << ",";
		outputFile << countData.movedAdjacent[x] << ",";
		outputFile << countData.movedConserved[x] << ",";
		outputFile << countData.mutualConserved[x] << endl;
	}
}





