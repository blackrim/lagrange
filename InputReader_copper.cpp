/*
 * InputReader.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#include <string>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iostream>
using namespace std;

#include "InputReader_copper.h"
#include "BioGeoTreeTools_copper.h"
#include "Utils.h"
#include "RateMatrixUtils.h"

#include "tree_reader.h"
#include "tree.h"
#include "node.h"

InputReader::InputReader():nareas(0),nspecies(0){
}

void InputReader::readMultipleTreeFile(string filename, vector<Tree *> & ret){
	TreeReader tr;
	ifstream ifs( filename.c_str() );
	string temp;
	int count = 1;
	while( getline( ifs, temp ) ){
		if(temp.size() > 1){
			Tree * intree = tr.readTree(temp);
			cout << "Tree "<< count <<" has " << intree->getExternalNodeCount() << " leaves." << endl;
			ret.push_back(intree);
			count++;
		}
	}
}

map<string,vector<int> > InputReader::readStandardInputData(string filename){
	ifstream ifs( filename.c_str() );
	string temp;
	bool first = false;
	nareas=0;
	nspecies=0;
	map<string,vector<int> > data;
	string line;
	vector<string> tokens;
	string del("\t ");
	while(getline(ifs,line)){
		if (first == false){
			first = true;
			tokens.clear();
			Tokenize(line, tokens, del);
			for(unsigned int j=0;j<tokens.size();j++){
				TrimSpaces(tokens[j]);
			}
			nspecies = atoi(tokens[0].c_str());
			nareas = atoi(tokens[1].c_str());
		}else{
			tokens.clear();
			Tokenize(line, tokens, del);
			for(unsigned int j=0;j<tokens.size();j++){
				TrimSpaces(tokens[j]);
			}
			cout << "Reading species: " << tokens[0] << " ";
			vector<int> speciesdata(nareas,0);
			for(int i=0;i<nareas;i++){
				char spot = tokens[1][i];
				if (spot == '1')
					speciesdata[i] = 1;
				cout << spot - '0';
			}
			cout << endl;
			data[tokens[0]] = speciesdata;
		}
	}
	ifs.close();
	return data;
}

void InputReader::checkData(map<string,vector<int> > data ,vector<Tree *> trees){
	vector<string> dataspecies;
	map<string,vector<int> >::const_iterator itr;
	for(itr = data.begin(); itr != data.end(); ++itr){
		dataspecies.push_back(itr->first);
	}
	//for(unsigned int i=0;i < trees.size();i++){
		for(int j=0;j<trees[0]->getExternalNodeCount();j++){
			bool test = false;
			for (unsigned int k=0;k<dataspecies.size();k++){
				if (trees[0]->getExternalNode(j)->getName() == dataspecies[k])
					test = true;
			}
			if(test == false){
				cout << "Error: " << trees[0]->getExternalNode(j)->getName() << " not found in data file." << endl;
				exit(0);
			}
		}
	//}
}
