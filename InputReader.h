/*
 * InputReader.h
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#ifndef INPUTREADER_COPPER_H_
#define INPUTREADER_COPPER_H_

#include <vector>
#include <string>
using namespace std;

#include "tree.h"

class InputReader{
	public:
	InputReader();
		void readMultipleTreeFile(string filename,vector<Tree *>&);
		map<string,vector<int> > readStandardInputData(string filename);
		void checkData(map<string,vector<int> >,vector<Tree *>);
		int nareas;
		int nspecies;
};

#endif /* INPUTREADER_H_ */
