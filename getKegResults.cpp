#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <utility>


///////////////////////////////////////////////////
///CONVERT ALL IDS TAGS TO UPPERCASE WHEN PARSING//
///////////////////////////////////////////////////

using namespace std;

struct geneInfo{ //stores parsed info from genbank
	string oldLocusTag;
	string locusTag;
	string proteinID;
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

void parseGenBank(ifstream& genBankFile, vector<geneInfo>& parsedInfo);
void parseKegFile(ifstream& kegFile, vector<kegInfo>& parsedKeg);
void parseSyntenyResults(ifstream& syntenyResults, vector<syntenyResult>& parsedInfo);
void removeMismatches(vector<syntenyResult>& forward, vector<syntenyResult>& reverse);
void findMutualConserved(vector<syntenyResult>& forward, vector<syntenyResult>& reverse);
void sortResults(vector<syntenyResult> forward, vector<syntenyResult> reverse, sortedResults& results);
void categorizeResults(sortedResults results, vector<geneInfo> parsedGenBank, vector<kegInfo> parsedKeg, ofstream& outputFile);
vector<int> getCategoryCounts(vector<pair<string, string> > results, vector<geneInfo> parsedGenBank, vector<kegInfo> parsedKeg, ofstream& outputFile);
string parseValue(string line);
string parseKegLine(string line);

int main(int argc, char *argv[]){
	if (argc != 7){
		cout << "missing/too many arguments!"<< endl;
		return 0;
	}
	
	ifstream genBankFile, kegFile ,forwardSyntenyFile, reverseSyntenyFile, kegLabelFile;
	ofstream countFile;
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
	
	kegLabelFile.open(argv[5]);
	if (!kegLabelFile.is_open()){
		ofstream kegLabelsFile;
		kegLabelsFile.open(argv[5]);
		for(int x = 0; x < kegParsed.size(); x++){
			kegLabelsFile << kegParsed[x].category<<",";
		}
		kegLabelsFile << endl;
	}else{
		cout << "KEG label file found" << endl;
	}
	
	countFile.open(argv[6]);
	categorizeResults(resultsSorted, genBankParsed, kegParsed, countFile);


	
	cout << "UNMOVED: " << resultsSorted.not_moved.size() <<endl;
	cout << "MOVED: " << resultsSorted.moved.size() <<endl;
	cout << "MOVED ALONE: " << resultsSorted.moved_adjacent.size() <<endl;
	cout << "MOVED INTO CONSERVED REGION: " << resultsSorted.moved_conserved.size() <<endl;
	cout << "MOVED INTO/FROM CONSERVED REGION: " << resultsSorted.conserved_both.size() <<endl;
	
	return 0;
}

//parses out the locus tag, old locus tag, and the corresponding protein id from the genbank
//this will be used to link the keg information and the synteny results
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

//finds categories and the genes within those categories
//assigns them to values of a struct then adds the struct to a vector
void parseKegFile(ifstream& kegFile, vector<kegInfo>& parsedKeg){
	string line;
	kegInfo keg;
	bool passedOverview = false;
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
			keg.genes.push_back(parseKegLine(line));
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
	
	return parsedLine;
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

//any mismatches between the forward and reverse results will be removed
//this does not filter out mismatches between movement into/from conserved regions
void removeMismatches(vector<syntenyResult>& forward, vector<syntenyResult>& reverse){
	int count = 0;
	int matchCount;
	vector<syntenyResult> forwardTemp, reverseTemp;
	//cout << "total forward: " << forward.size() <<endl;
	//cout << "total reverse: " << reverse.size() <<endl;
	
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
	
	cout << "remaining forward: " << forward.size() << endl;
	cout << "remaining reverse: " << reverse.size() << endl;
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
					//cout <<forward[i].sub_prot << " " << reverse[j].sub_prot <<endl;
					count++;
				}
			}
		}
	}
	//cout << "mutually moved into/from conserved: " << count <<endl;
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
		/*if (forward[i].moved ==1 && forward[i].moved_adjacent == 1 && forward[i].moved_conserved == 0){ //moved based upon adjacent proteins
			temp.first = forward[i].sub_prot;
			temp.second = forward[i].query_prot;
			results.moved_adjacent.push_back(temp);
		}*/
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
	vector<int> notMoved, movedAbsolute, movedAdjacent, movedConserved, moveConservedBoth;
	notMoved = getCategoryCounts(results.not_moved, parsedGenBank, parsedKeg, outputFile);
	movedAbsolute = getCategoryCounts(results.moved, parsedGenBank, parsedKeg, outputFile);
	movedAdjacent = getCategoryCounts(results.moved_adjacent, parsedGenBank, parsedKeg, outputFile);
	movedConserved = getCategoryCounts(results.moved_conserved, parsedGenBank, parsedKeg, outputFile);
	movedAdjacent = getCategoryCounts(results.conserved_both, parsedGenBank, parsedKeg, outputFile);
	outputFile << endl;
}


//I'm writing this to be compatable with locus tags or old locus tags
vector<int> getCategoryCounts(vector<pair<string, string> > results, vector<geneInfo> parsedGenBank, vector<kegInfo> parsedKeg, ofstream& outputFile){
	vector<pair<string, string> > locusTags;
	vector<int> categoryCounts;
	pair<string, string> temp;
	for(int x = 0; x < results.size(); x++){
		for (int y=0; y < parsedGenBank.size(); y ++){
			if(results[x].first == parsedGenBank[y].proteinID){
				temp.first = "";
				temp.second = "";
				//cout << results[x].first << "\t" << parsedGenBank[y].oldLocusTag << "\t" << parsedGenBank[y].locusTag <<endl;
				temp.first = parsedGenBank[y].oldLocusTag;
				temp.second = parsedGenBank[y].locusTag;
				locusTags.push_back(temp);
			}
		}
	}
	
	int count;
	for(int i =0; i < parsedKeg.size(); i++){
		count = 0;
		for(int j = 0; j < parsedKeg[i].genes.size(); j++){
			for(int k = 0; k < locusTags.size(); k++){
				if(parsedKeg[i].genes[j] == locusTags[k].first || parsedKeg[i].genes[j] == locusTags[k].second){
					count++;
				}
			}
		}
		categoryCounts.push_back(count);
		//cout << parsedKeg[i].category << "\t" << count <<endl;
		outputFile << count << ",";
	}
	return categoryCounts;
}







