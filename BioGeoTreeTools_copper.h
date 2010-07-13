/*
 * PhyloTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

#ifndef PHYLOTREE_COPPER_H_
#define PHYLOTREE_COPPER_H_

#include <string>
#include <map>
#include <vector>
#include "AncSplit.h"
using namespace std;

#include "tree.h"
#include "node.h"

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

class BioGeoTreeTools_copper {
public :
	Tree * getTreeFromString(string treestring);
	vector<Node *> getAncestors(Tree & tree, Node & node);

	void summarizeSplits(Node * node,map<vector<int>,vector<AncSplit> > & ans,map<int,string> &areanamemaprev, RateModel * rm);
#ifdef BIGTREE
	void summarizeAncState(Node * node,vector<mpfr_class> & ans,map<int,string> &areanamemaprev, RateModel * rm);
#else
	void summarizeAncState(Node * node,vector<double> & ans,map<int,string> &areanamemaprev, RateModel * rm);
#endif
	string get_string_from_dist_int(int dist,map<int,string> &areanamemaprev, RateModel * rm);
};

#endif /* PHYLOTREE_H_ */
