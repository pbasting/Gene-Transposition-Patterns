/***************************************************************************************************
FormatKegResults
Written by: Preston Basting
Email:pjb68507@uga.edu
Lab: Jan Mrazek
Last Changed: 9/18/2017
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
	vector<int> notMoved, moved, movedConserved, mutualConserved;
};

struct movements{
	string subject;
	string query;
	string keg;
	int move;
};


//takes a keg file and parses out the categories. Stores all the categories into a vector
vector<string> parseKegFile(ifstream& kegFile);

//takes a line from a keg file that contains a category. Parses out the name of the category
string parseKegLine(string line);

//converts a string to uppercase
string upperCase(string line);


//takes the concatenated genus results from 'getKegResults', parses and stores in a struct
void buildMovementResults(ifstream& data, vector<movements>& proteins);


//takes a parsed line from the concatenated genus results
//returns an integer corresponding to the movement category
int getMovementCategory(string line);

//called by buildMovementResults
//takes a line from the concatenated genus results
//parses the line into a struct
void buildResult(string line, movements& result);

//takes the vector of structs containing all the parsed data for a genus
//erases the values of any duplicates found
//keeps the pair with the highest movement classification
void removeDuplicates(vector<movements>& proteins);

//takes the parsed results from 'getKegResults', counts how many proteins match to each category and splits
//the count into movement categories
void buildTable(vector<string> categories, vector<movements> results, data& countData);

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
	
	vector<movements> movementResults;
	buildMovementResults(kegResults, movementResults);
	removeDuplicates(movementResults);
	
	kegResults.clear();
	kegResults.seekg(0,ios::beg);
	data countData;
	buildTable(kegCategories, movementResults, countData);
	
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

void buildMovementResults(ifstream& data, vector<movements>& proteins){
	string line = "";
	string temp="";
	movements result;
	int pos;
	while(!data.eof()){
		getline(data, line);
		if(line.find("!!")!= string::npos){ //indicates movement category line
			result.move = getMovementCategory(line);
			while(line.find("**")==string::npos){ //'**' indicates the end of the results
				getline(data, line);
				if (line.find("$$")!=string::npos){ //'$$' indicates a line containing data
					buildResult(line, result);
					proteins.push_back(result);
				}
			}
		}
	}
}


int getMovementCategory(string line){
	if (line.find("NOT_MOVED")!=string::npos){
		return 0;
	}
	if (line.find("MOVED_ADJACENT")!=string::npos){
		return 1;
	}
	if (line.find("MOVED_CONSERVED")!=string::npos){
		return 2;
	}
	if (line.find("MOVED_MUTUAL_CONSERVED")!=string::npos){
		return 3;
	}
}

void buildResult(string line, movements& result){
	int pos=0;
	int endPos;
	string temp = "";
	line =  line.substr(line.find("\t")+1); //moves past first tab
	endPos = line.find("\t");
	while(pos < endPos){
		temp+=line[pos];
		pos++;
	}
	result.subject = temp;
	line = line.substr(line.find("\t")+1);//moves past subject
	pos=0;
	endPos = line.find("\t");
	temp = "";
	while(pos < endPos){
		temp+=line[pos];
		pos++;
	}
	result.query = temp;
	line = line.substr(line.find("\t")+1); //moves past query
	line = line.substr(0,line.find("/product=")); //doesnt include genbank product info
	result.keg = line;
}


void removeDuplicates(vector<movements>& proteins){
	int count = 0;
	int total = 0;
	int y;
	cout << "REMOVING DUPLICATES..." <<endl;
	for(int x=0; x < proteins.size(); x++){ //loops through all proteins
		total++;
		if(proteins[x].subject.length() > 0){ //doesnt search for matches if already erased
			y = x+1;
			while(y < proteins.size()){ //loops through all downstream proteins //upstream ones have already been checked
				//looks for a match
				if ((proteins[x].subject == proteins[y].subject || proteins[x].subject == proteins[y].query || 
				proteins[x].query == proteins[y].subject || proteins[x].query == proteins[y].subject)&& (x != y)
				&&(proteins[x].subject.length() > 0) && (proteins[x].query.length() > 0)){
					//erases the pair with the lower move category
					if (proteins[y].move > proteins[x].move){
							proteins[x].subject = "";
							proteins[x].query = "";
							proteins[x].move = -1;
							proteins[x].keg = "";
							count++;
							if (count%1000 == 0){
								cout << count << " removed" <<endl;
							}
							break;
					}else{
							proteins[y].subject = "";
							proteins[y].query = "";
							proteins[y].move = -1;
							proteins[y].keg = "";
							count++;
							if (count%1000 == 0){
								cout << count << " removed" <<endl;
							}
					}
				}
				y++;
			}
		
		}
	
	}
	cout << count << "out of " << total << " total Protein pairs" <<endl;
}

void buildTable(vector<string> categories, vector<movements> results, data& countData){
	countData.categories = categories;
	//adds zeros to all vector positions so specific indexes can be increased during the count
	for(int i = 0; i < categories.size(); i++){
			countData.notMoved.push_back(0);
			countData.moved.push_back(0);
			countData.movedConserved.push_back(0);
			countData.mutualConserved.push_back(0);
	}
	for(int x = 0; x < results.size(); x++){
		for(int y = 0; y < categories.size(); y++){
			//tallies matches at the indices that correspond with the category
			if (results[x].keg.find(categories[y])!=string::npos){
				switch (results[x].move){
					case-1 : break;
					case 0 : countData.notMoved[y]++;
							 break;
					case 1 : countData.moved[y]++;
							 break;
					case 2 : countData.movedConserved[y]++;
							 break;
					case 3 : countData.mutualConserved[y]++;
							 break;
				
				}
			}
		}
	}
}


void outputTable(data countData, ofstream& outputFile){
	outputFile << "FUNCTION,UNMOVED,MOVED,MOVED.CONS,MUTUAL.CONS"<<endl;
	for (int x = 0; x < countData.categories.size(); x++){
		if(countData.categories[x] == upperCase("Folding, sorting and degradation")){ //the comma in this category messes up the comma delimiting
			outputFile << upperCase("Folding sorting and degradation,"); //this is the same catagory title without the comma
		}else{
			outputFile << countData.categories[x] << ",";
		}
		outputFile << countData.notMoved[x] << ",";
		outputFile << countData.moved[x] << ",";
		outputFile << countData.movedConserved[x] << ",";
		outputFile << countData.mutualConserved[x] << endl;
	}
}





