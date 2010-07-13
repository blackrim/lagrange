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
using namespace std;

#include "InputReader.h"
#include "BioGeoTreeTools.h"
#include "Utils.h"
#include "RateMatrixUtils.h"

#include <Phyl/TreeTemplate.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Utils/BppVector.h>
using namespace bpp;

InputReader::InputReader(){
}

vector<TreeTemplate<Node> *> InputReader::readMultipleTreeFile(string filename){
	BioGeoTreeTools * nr = new BioGeoTreeTools();
	ifstream ifs( filename.c_str() );
	string temp;
	vector<TreeTemplate<Node> *> ret;
	int count = 1;
	while( getline( ifs, temp ) ){
		if(temp.size() > 1){
			TreeTemplate<Node> * intree = nr->getTreeFromString(temp);
			cout << "Tree "<< count <<" has " << intree->getNumberOfLeaves() << " leaves." << endl;
			ret.push_back(intree);
			count++;
		}
	}
	delete nr;
	return ret;
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

void InputReader::checkData(map<string,vector<int> > data,vector<TreeTemplate<Node> *> trees){
	vector<string> dataspecies;
	map<string,vector<int> >::const_iterator itr;
	for(itr = data.begin(); itr != data.end(); ++itr){
		dataspecies.push_back(itr->first);
	}
	for(unsigned int i=0;i < trees.size();i++){
		vector<string> lnames = trees[0]->getLeavesNames();
		for(unsigned int j=0;j<lnames.size();j++){
			bool test = false;
			for (unsigned int k=0;k<dataspecies.size();k++){
				if (lnames[j] == dataspecies[k])
					test = true;
			}
			if(test == false){
				cout << "Error: " << lnames[j] << " not found in data file." << endl;
				exit(0);
			}
		}
	}
}
