/*
 * PhyloTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

#ifndef PHYLOTREE_H_
#define PHYLOTREE_H_

#include <string>
#include <map>
#include <vector>
#include <Phyl/Node.h>
#include <Phyl/TreeTemplate.h>
#include "AncSplit.h"
using namespace bpp;
using namespace std;

class BioGeoTreeTools {
	public :
		TreeTemplate<Node> * getTreeFromString(string treestring) throw (Exception);
		vector<int> getAncestors(TreeTemplate<bpp::Node> & tree, int nodeId);
		int getLastCommonAncestor(TreeTemplate<bpp::Node> & tree, const vector<int>& nodeIds);
		void summarizeSplits(Node * node,map<vector<int>,vector<AncSplit> > & ans,map<int,string> &areanamemaprev, RateModel * rm);
		void summarizeAncState(Node * node,vector<double> & ans,map<int,string> &areanamemaprev, RateModel * rm);
};

#endif /* PHYLOTREE_H_ */
