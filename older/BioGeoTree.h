/*
 * BioGeoTree.h
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */

#ifndef BIOGEOTREE_H_
#define BIOGEOTREE_H_

#include <vector>
#include <string>
#include <map>
using namespace std;

#include "RateModel.h"
#include "AncSplit.h"
#include "BranchSegment.h"

#include <Phyl/TreeTemplateTools.h>
#include <Phyl/TreeTemplate.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Utils/BppVector.h>

class BioGeoTree{
private:
	bpp::TreeTemplate<bpp::Node> * tree;
	int numofnodes;
	vector<double> periods;
	string seg;
	string age;
	string dc;
	string en;
	string andc;
	//reverse bits
	string revB;
	//end reverse bits
	vector<int> * columns;
	vector<int> * whichcolumns;
	RateModel * rootratemodel;
	map<int, vector<int> > * distmap; // a map of int and dist
	bool store_p_matrices;
	bool use_stored_matrices;
	double scale;

	/*
	 * replaces tree.getNodeID
	 * when given i it will return tree->getNode(i)
	 */
	map<int,bpp::Node *> tree_get_node_from_id;
	/*
	 * benchmark variables
	 */
	clock_t cl1;
	clock_t cl2;
	clock_t c3;
	clock_t c4;
	clock_t c5;
	clock_t c6;

public:
	BioGeoTree(bpp::TreeTemplate<bpp::Node> * tr, vector<double> ps);
	void set_store_p_matrices(bool);
	void set_use_stored_matrices(bool);
	void set_default_model(RateModel * mod);
	void update_default_model(RateModel * mod);
	double eval_likelihood(bool marg);
	void set_excluded_dist(vector<int> ind,bpp::Node * node);
	void set_tip_conditionals(map<string,vector<int> > distrib_data);
	bpp::Vector<double> conditionals(bpp::Node & node, bool marg, bool sparse);
	void ancdist_conditional_lh(bpp::Node & node, bool marg);

/*
	fossil data
 */
	void setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea);
	void setFossilatNodeByMRCA_id(int id, int fossilarea);
	void setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea, double age);
	void setFossilatBranchByMRCA_id(int id, int fossilarea, double age);

/*
	for calculating forward and reverse
 */
	void prepare_ancstate_reverse();
	void reverse(bpp::Node &);
	map<vector<int>,vector<AncSplit> > calculate_ancsplit_reverse(bpp::Node & node,bool marg);
	vector<double> calculate_ancstate_reverse(bpp::Node & node,bool marg);

/*
	for timing things
 */
	double ti;
	double ti2;
	double ti3;
};

#endif /* BIOGEOTREE_H_ */
