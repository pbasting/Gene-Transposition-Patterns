/***************************************************************************************************
getKegResults
Written by: Preston Basting
Email:pjb68507@uga.edu
Lab: Jan Mrazek
Last Changed: 9/1/2017
Purpose: This is a component of a series of programs designed to classify protein
		 'movement' when comparing two organisms and determine if proteins belonging
		 to different functional categories are more likely to 'move'
		 
		 This function takes the results from 'CompareOrthologs', removes all missmatches from the 
		 forward and reverse results, converts the protein IDs to locus tags found in the genbank file,
		 then classifies the proteins into functional categories based upon the .kegg file from genome.jp.
		 The results are output to a text file for future concatenation with results from other organisms
		 within the same genus.
		 
Arguments: (1)subject genbank file, (2)subject keg file, (3)forward 'CompareOrthologs' results
		   (4)reverse 'CompareOrthologs' results, (5) name of output file
****************************************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <utility>
#include <algorithm>


using namespace std;

struct geneInfo{ //stores parsed info from genbank
	string oldLocusTag;
	string locusTag;
	string proteinID;
	string product;
};

struct syntenyResult{ //stores parsed info from synteny results
	string sub_prot;
	string query_prot;
	int moved;
	int moved_adjacent;
	int moved_conserved;
	int conserved_both;
};

struct sortedResults{ //stores synteny results sorted by movement category
	vector<pair<string, string> > not_moved;
	vector<pair<string, string> > moved;
	vector<pair<string, string> > moved_adjacent;
	vector<pair<string, string> > moved_conserved;
	vector<pair<string, string> > conserved_both;
};

struct kegInfo{ //stores parsed info from keg
	string category;
	vector<string> genes;
};

//takes subject genbank file and parses the relevant information into a vector of structs
void parseGenBank(ifstream& genBankFile, vector<geneInfo>& parsedInfo);

//takes the subject keg file and parses the relavent information into a vector of structs
void parseKegFile(ifstream& kegFile, vector<kegInfo>& parsedKeg);

//takes the results from 'CompareOrthologs' and parses the information into a vector of structs
void parseSyntenyResults(ifstream& syntenyResults, vector<syntenyResult>& parsedInfo);

//takes the parsed information from the forward and reverse synteny results and removes any
//mismatches or matches in one result that don't exist in the other
void removeMismatches(vector<syntenyResult>& forward, vector<syntenyResult>& reverse);

//takes the paresd forward and reverse synteny results and finds the proteins that enter
//conserved regions from the forward and reverse perspective
void findMutualConserved(vector<syntenyResult>& forward, vector<syntenyResult>& reverse);

//sorts the results into the movement categories and stores that as a struct
void sortResults(vector<syntenyResult> forward, vector<syntenyResult> reverse, sortedResults& results);

//assigns keg functional categories to the parsed results and outputs to a text file
//uses the genbank to convert protein IDs to locus tags which are used by the keg file
void categorizeResults(sortedResults results, vector<geneInfo> parsedGenBank, vector<kegInfo> parsedKeg, ofstream& outputFile);

//used by 'categorizeResults' function
//finds the appropriate keg category (if it exists) and exports the subject and query protein IDs
//and the keg categories to a text file
void getCategoryCounts(vector<pair<string, string> > results, vector<geneInfo> parsedGenBank, vector<kegInfo> parsedKeg, ofstream& outputFile);

string parseValue(string line);
string parseKegLine(string line);
string upperCase(string line);
string removePosition(string proteinID);

int main(int argc, char *argv[]){
	if (argc != 6){
		cout << "missing/too many arguments!"<< endl;
		return 0;
	}
	
	ifstream genBankFile, kegFile ,forwardSyntenyFile, reverseSyntenyFile, kegLabelFile;
	ofstream countFile;
	genBankFile.open(argv[1]);
	kegFile.open(argv[2]);
	forwardSyntenyFile.open(argv[3]);
	reverseSyntenyFile.open(argv[4]);
	countFile.open(argv[5]);
	
	vector<geneInfo> genBankParsed;
	vector<kegInfo> kegParsed;
	vector<syntenyResult> forwardResults, reverseResults;
	sortedResults resultsSorted;
	
	//parses input files into structs
	parseGenBank(genBankFile, genBankParsed);
	parseKegFile(kegFile, kegParsed);
	parseSyntenyResults(forwardSyntenyFile, forwardResults);
	parseSyntenyResults(reverseSyntenyFile, reverseResults);
	
	//cleans up the synteny data
	removeMismatches(forwardResults, reverseResults);
	findMutualConserved(forwardResults, reverseResults);
	sortResults(forwardResults, reverseResults, resultsSorted);
	
	//gets title
	string title = argv[3];
	title = title.substr((title.rfind("/")+1)); //removes path from title
	
	//outputs count information to a file
	countFile <<endl<<endl<<endl<<endl<< "/////////////////////////////////////////////////" << endl <<endl;
	countFile << "##" << upperCase(title) << endl;
	countFile <<endl<< "/////////////////////////////////////////////////"<<endl;
	countFile << "TOTAL: " << forwardResults.size() <<endl;
	countFile << "NOT MOVED: " << resultsSorted.not_moved.size() <<endl;
	countFile << "MOVED ABSOLUTE: " << resultsSorted.moved.size() <<endl;
	countFile<< "MOVED ADJACENT: " << resultsSorted.moved_adjacent.size() <<endl;
	countFile << "MOVED CONSERVED: " << resultsSorted.moved_conserved.size() <<endl;
	countFile << "MOVED MUTUAL CONSERVED: " << resultsSorted.conserved_both.size() <<endl;
	
	//outputs protein IDs and keg categories to a file
	categorizeResults(resultsSorted, genBankParsed, kegParsed, countFile);
/*
	cout << "UNMOVED: " << resultsSorted.not_moved.size() <<endl;
	cout << "MOVED: " << resultsSorted.moved.size() <<endl;
	cout << "MOVED ALONE: " << resultsSorted.moved_adjacent.size() <<endl;
	cout << "MOVED INTO CONSERVED REGION: " << resultsSorted.moved_conserved.size() <<endl;
	cout << "MOVED INTO/FROM CONSERVED REGION: " << resultsSorted.conserved_both.size() <<endl;
*/
	return 0;
}

//parses out the locus tag, old locus tag, and the corresponding protein id from the genbank
//this will be used to link the keg information and the synteny results
void parseGenBank(ifstream& genBankFile, vector<geneInfo>& parsedInfo){
	string line;
	geneInfo gene;
	while(!genBankFile.eof()){
		getline(genBankFile, line);
		if (line.find("    CDS    ")!=string::npos){ //finds CDS
			gene.locusTag = "";
			gene.oldLocusTag = "";
			gene.proteinID="";
			gene.product="";
			do{
				getline(genBankFile, line);
				if(line.find("/locus_tag=") != string::npos){ //finds locus tag
					gene.locusTag = parseValue(line);
				}
				if(line.find("/old_locus_tag=") != string::npos){ //finds locus tag
					gene.oldLocusTag = parseValue(line);
				}
				if(line.find("/product=") != string::npos){
					gene.product = ("/product="+parseValue(line));
				}
				if(line.find("/protein_id=") != string::npos){ //find protein ID
					gene.proteinID = parseValue(line);
					gene.proteinID = gene.proteinID.substr(0, gene.proteinID.rfind(".")); //removes version number
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
	for (pos; ((line[pos] != '\n') && (line[pos] != '"')) && (value.find("\" ") ==string::npos) && pos < line.length(); pos++){
		value+=line[pos];
	}
	return upperCase(value);
}

//finds categories and the genes within those categories
//assigns them to values of a struct then adds the struct to a vector
void parseKegFile(ifstream& kegFile, vector<kegInfo>& parsedKeg){
	string line, temp;
	kegInfo keg;
	bool passedOverview = false;
	bool match = false;
	while(!kegFile.eof()){
		getline(kegFile, line);
		if(line.find("<b>Overview</b>")==string::npos && line[0] == 'B' && line.find("<b>") != string::npos){
			if(passedOverview){ //first category is overview which is skipped
				parsedKeg.push_back(keg);
				keg.category = "";
				keg.genes.clear();
			}
			keg.category = parseKegLine(line);
			passedOverview = true; //assigned to true after overview is skipped
		}
		if(line[0] == 'D'){
			match = false;
			temp = parseKegLine(line);
			for(int i = 0; i < keg.genes.size(); i++){
				if (keg.genes[i] == temp){
					match = true;
				}
			}
			if (!match){//prevents duplicates from being added
				keg.genes.push_back(temp); 
			}
		}
	}	
}

string parseKegLine(string line){
	int pos;
	string parsedLine = "";
	if(line.find("<b>") != string::npos){ //if the line provided is a category
		pos = line.find("<b>");
		pos+=3;
		while(line[pos] != '<' && pos < line.length()){
			parsedLine += line[pos];
			pos++;
		}
	}else{ //if the line provided is a gene
		pos = 7;
		while (line[pos] != ' ' && pos < line.length()){
			parsedLine += line[pos];
			pos++;
		}
	}
	
	return upperCase(parsedLine);
}

//parses the results from 'CompareOrthologs' and stores them into a struct
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
			for (pos; ((line[pos] != '\n') && (line[pos-1] != ',')); pos++){
				temp += line[pos]; 
			}
			temp = temp.substr(0,temp.rfind(","));
			result.sub_prot = upperCase(temp);
			
			//gets query protein ID///
			temp = "";
			pos = line.rfind("_prot_");
			pos+=6;
			for (pos; ((line[pos] != '\n') && (line[pos-1] != ',')); pos++){
				temp += line[pos];
			}
			temp = temp.substr(0,temp.rfind(","));
			result.query_prot = upperCase(temp);
			
			//gets movement info//
			line = line.substr(line.find(temp));
			line = line.substr(line.find(",")+1); //moves past query protein
			line = line.substr(line.find(",")+1); // moves past subject protein pos
			line = line.substr(line.find(",")+1); // moves past query protein pos
			result.moved = line[0]-48; //-48 to convert char to int
			line = line.substr(line.find(",")+1); //moves past 'moved'
			result.moved_adjacent = line[0]-48; //-48 to convert char to int
			line = line.substr(line.find(",")+1); //moves past 'moved adjacent'
			result.moved_conserved = line[0]-48; //-48 to convert char to int
			result.conserved_both = 0;
			parsedInfo.push_back(result);
		}
	}
}

//any mismatches between the forward and reverse results will be removed
//this does not filter out mismatches between movement into/from conserved regions
void removeMismatches(vector<syntenyResult>& forward, vector<syntenyResult>& reverse){
	int count = 0;
	int matchCount;
	vector<syntenyResult> forwardTemp, reverseTemp;
	
	//remove proteins with no matches in the reverse condition//
	if(forward.size() < reverse.size()){
		for (int x=0; x < forward.size(); x++){
			for (int y=0; y < reverse.size(); y++){
				if (forward[x].sub_prot == reverse[y].query_prot){
					forwardTemp.push_back(forward[x]);
					reverseTemp.push_back(reverse[y]);
				}
			}
		}
	}else{
		for (int x=0; x < reverse.size(); x++){
			for (int y=0; y < forward.size(); y++){
				if (reverse[x].sub_prot == forward[y].query_prot){
					reverseTemp.push_back(reverse[x]);
					forwardTemp.push_back(forward[y]);
				}
			}
		}
	}
	forward.clear();
	reverse.clear();
	forward = forwardTemp;
	reverse = reverseTemp;
	forwardTemp.clear();
	reverseTemp.clear();
	
	//removes mismatches and and pairs with divergence in movement classification
	
	for(int i = 0; i < forward.size(); i++){
		for (int j =0; j < reverse.size(); j++){
			if (forward[i].sub_prot == reverse[j].query_prot && forward[i].query_prot == reverse[j].sub_prot &&
				forward[i].moved == reverse[j].moved && forward[i].moved_adjacent == reverse[j].moved_adjacent){
				forwardTemp.push_back(forward[i]);
				reverseTemp.push_back(reverse[j]);
			}
		}
	}
	
	forward.clear();
	reverse.clear();
	forward = forwardTemp;
	reverse = reverseTemp;
}


//determines if the forward and reverse results match regarding the movement into a conserved region
// if they match the variable 'conserved_both' is set to 1
void findMutualConserved(vector<syntenyResult>& forward, vector<syntenyResult>& reverse){
	int count = 0;
	for (int i =0; i < forward.size(); i++){
		for (int j = 0; j < reverse.size(); j++){
			if (forward[i].sub_prot == reverse[j].query_prot){
				if (forward[i].moved ==1 && forward[i].moved_adjacent ==1 &&forward[i].moved_conserved ==1 &&
				    reverse[j].moved ==1 && reverse[j].moved_adjacent ==1 &&reverse[j].moved_conserved ==1){
					
					forward[i].conserved_both = 1;
					reverse[j].conserved_both = 1;
					count++;
				}
			}
		}
	}
}


//sorts the synteny results into 5 categories and stores them in a struct
void sortResults(vector<syntenyResult> forward, vector<syntenyResult> reverse, sortedResults& results){
	pair<string, string> temp;
	for (int i = 0; i < forward.size(); i++){
		if (forward[i].moved == 0){ //didnt move
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.not_moved.push_back(temp);
		}
		if (forward[i].moved ==1 && forward[i].moved_adjacent == 0){ //moved based on overall position but not adjacent proteins
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.moved.push_back(temp);
		}
		if (forward[i].moved ==1 && forward[i].moved_adjacent == 1 && forward[i].moved_conserved == 1 && forward[i].conserved_both == 0){ //moved based upon adjacent proteins into a conserved region
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.moved_conserved.push_back(temp);
		}
		//moved based upon adjacent proteins from a conserved region into a conserved region
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
		
		//gets pairs that have moved based upon adjacent proteins but haven't entered or left a conserved zone
	}
	bool match = false;
	for (int i = 0; i < forward.size(); i++){
		if (forward[i].moved ==1 && forward[i].moved_adjacent == 1 && forward[i].moved_conserved == 0){
			match = false;
			for(int x = 0; x < results.moved_conserved.size(); x++){
				if (results.moved_conserved[x].first == forward[i].sub_prot){
					match = true;
				}
			}
			if (!match){
				temp.first = forward[i].sub_prot;
				temp.second = forward[i].query_prot;
				results.moved_adjacent.push_back(temp);
			}
		}
	}
}


void categorizeResults(sortedResults results, vector<geneInfo> parsedGenBank, vector<kegInfo> parsedKeg, ofstream& outputFile){
	outputFile << "!!NOT_MOVED!!"<<endl;
	getCategoryCounts(results.not_moved, parsedGenBank, parsedKeg, outputFile);
	outputFile << "**" <<endl;
	outputFile << "!!MOVED_ABSOLUTE!!"<<endl;
	getCategoryCounts(results.moved, parsedGenBank, parsedKeg, outputFile);
	outputFile << "**" <<endl;
	outputFile << "!!MOVED_ADJACENT!!"<<endl;
	getCategoryCounts(results.moved_adjacent, parsedGenBank, parsedKeg, outputFile);
	outputFile << "**" <<endl;
	outputFile << "!!MOVED_CONSERVED!!"<<endl;
	getCategoryCounts(results.moved_conserved, parsedGenBank, parsedKeg, outputFile);
	outputFile << "**" <<endl;
	outputFile << "!!MOVED_MUTUAL_CONSERVED!!"<<endl;
	getCategoryCounts(results.conserved_both, parsedGenBank, parsedKeg, outputFile);
	outputFile << "**" <<endl;
}


void getCategoryCounts(vector<pair<string, string> > results, vector<geneInfo> parsedGenBank, vector<kegInfo> parsedKeg, ofstream& outputFile){
	bool categorized = false;
	int productIndex=0;
	bool productFound = false;
	for(int x=0; x< results.size(); x++){ //loops through protein pairs
		outputFile <<"$$\t"<< results[x].first <<"\t" <<results[x].second << "\t";
		productIndex = 0;
		productFound = false;
		for(int y=0; y < parsedGenBank.size(); y++){
			if(removePosition(results[x].first) == parsedGenBank[y].proteinID){ //finds locus tag
				productIndex = y;
				productFound = true;
				for(int z=0; z < parsedKeg.size(); z++){ //loops through all keg categories
					for(int a=0; a < parsedKeg[z].genes.size(); a++){ //loops though all locus tags for a category
						
						if(parsedKeg[z].genes[a] == parsedGenBank[y].oldLocusTag || parsedKeg[z].genes[a] == parsedGenBank[y].locusTag){
							categorized = true; //controls whether or not a protein is assigned 'UNCATEGORIZED'
							outputFile << parsedKeg[z].category << "\t";
						}
					}
				}
				break;
			}
		}
		if(!categorized){
			outputFile << "UNCATEGORIZED" << "\t";
		}
		if(productFound == true){
			outputFile<< parsedGenBank[productIndex].product <<endl;
		}else{
			outputFile <<endl;
		}
		
		//outputFile <<endl;
		categorized = false;
	}
}

string upperCase(string line){
	transform(line.begin(), line.end(), line.begin(), ::toupper);
	return line;
}

string removePosition(string proteinID){
	proteinID = proteinID.substr(0,proteinID.find("."));
	return proteinID;
}






