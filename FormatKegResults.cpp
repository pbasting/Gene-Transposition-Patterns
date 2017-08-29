/***************************************************************************************************
FormatKegResults
Written by: Preston Basting
Lab: Jan Mrazek
Last Changed: 8/29/2017
Purpose: This is a component of a series of programs designed to classify protein
		 'movement' when comparing two organisms and determine if proteins belonging
		 to different functional categories are more likely to 'move'
		 
		 This function takes the concatenated 'getKegResults' results of all of the comparisions
		 in a genus and counts the number of hits for each kegg protein category based upon the movement 
		 category generates a table where each row is a different kegg category and each column is a 
		 different movement category. The values are the number or proteins that meat both classifications.
		 
Arguments: (1)any keg file, (2)concatenated genus results from 'getKegResults'
****************************************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <utility>
#include <algorithm>


using namespace std;

vector<string> parseKegFile(ifstream& kegFile);
string parseKegLine(string line);
string upperCase(string line);
vector<vector<string> > buildTable(vector<string> categories, ifstream& results);
////////////////////////////////////////////////////////////////////////////////
/////////////////develop notations for easy parsing in the keg results output///
////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){
	if (argc != 3){
		cout << "missing/too many arguments!"<< endl;
		return 0;
	}
	ifstream kegFile;
	kegFile.open(argv[1]);
	
	vector<vector<string> > data;
	vector<string> kegCategories = parseKegFile(kegFile);
	data.push_back(kegCategories);

	
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
		parsedLine += line[pos];
		pos++;
	}
	
	return upperCase(parsedLine);
}

string upperCase(string line){
	transform(line.begin(), line.end(), line.begin(), ::toupper);
	return line;
}

vector<vector<string> > buildTable(vector<string> categories, ifstream& results){
	///////////////////////////////////////////////////////
	/////////FINISH WRITING////////////////////////////////
	///////////////////////////////////////////////////////
}

