/*
 * InputReader.h
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#ifndef INPUTREADER_H_
#define INPUTREADER_H_

#include <vector>
#include <string>
using namespace std;

#include <Phyl/TreeTemplate.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Utils/BppVector.h>
using namespace bpp;

class InputReader{
	public:
		InputReader();
		vector<TreeTemplate<Node> *> readMultipleTreeFile(string filename);
		map<string,vector<int> > readStandardInputData(string filename);
		void checkData(map<string,vector<int> >,vector<TreeTemplate<Node> *>);
		int nareas;
		int nspecies;
};

#endif /* INPUTREADER_H_ */
