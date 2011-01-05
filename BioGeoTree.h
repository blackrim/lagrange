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

#include "tree.h"
#include "node.h"
#include "vector_node_object.h"

#include <armadillo>
using namespace arma;

//octave usage
//#include <octave/oct.h>

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

class BioGeoTree{
private:
	Tree * tree;
	vector<double> periods;
	string seg;
	string age;
	string dc;
	string en;
	string andc;
	vector<int> * columns;
	vector<int> * whichcolumns;
	RateModel * rootratemodel;
	map<int, vector<int> > * distmap; // a map of int and dist
	bool store_p_matrices;
	bool use_stored_matrices;

	//reverse bits
	string revB;
	bool rev;
	//end reverse bits

	//stochastic mapping bits
	string rev_exp_number;
	string rev_exp_time;
	bool stochastic;
	//map of period int and then branch length double
	map<int,map<double, mat > > stored_EN_matrices;
	map<int,map<double, mat > > stored_ER_matrices;
	//end mapping bits

	//scaling for larger trees
	double scale;
	double totalscale;
	double limit;
	bool run_with_scale;
	void scale_node(VectorNodeObject<double> * conds);

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
	BioGeoTree(Tree * tr, vector<double> ps);
	void set_store_p_matrices(bool);
	void set_use_stored_matrices(bool);
	void set_default_model(RateModel * mod);
	void update_default_model(RateModel * mod);
	double eval_likelihood(bool marg);
	void set_excluded_dist(vector<int> ind,Node * node);
	void set_tip_conditionals(map<string,vector<int> > distrib_data);
#ifdef BIGTREE
	VectorNodeObject<mpfr_class> conditionals(Node & node, bool marg, bool sparse);
#else
	VectorNodeObject<double> conditionals(Node & node, bool marg, bool sparse);
#endif
	//void ancdist_conditional_lh(bpp::Node & node, bool marg);
	void ancdist_conditional_lh(Node & node, bool marg);

/*
	fossil data
 */
	void setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea);
	void setFossilatNodeByMRCA_id(Node * id, int fossilarea);
	void setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea, double age);
	void setFossilatBranchByMRCA_id(Node * id, int fossilarea, double age);

/*
	for calculating forward and reverse
 */
	void prepare_ancstate_reverse();
	void reverse(Node &);
	map<vector<int>,vector<AncSplit> > calculate_ancsplit_reverse(Node & node,bool marg);
#ifdef BIGTREE
	vector<mpfr_class> calculate_ancstate_reverse(Node & node,bool marg);
#else
	vector<double> calculate_ancstate_reverse(Node & node,bool marg);
#endif
	~BioGeoTree();
	//need to override these at some point
	BioGeoTree(const BioGeoTree &L);             // copy constructor
    BioGeoTree & operator=(const BioGeoTree &L);

/*
 * for calculating forward and reverse for expected values (stochastic mapping)
 */
	void prepare_stochmap_reverse_all_nodes(int, int);
#ifdef BIGTREE
	vector<mpfr_class> calculate_reverse_stochmap(Node &, bool);
#else
	vector<double> calculate_reverse_stochmap(Node &, bool);
#endif
	vector<double> calculate_reverse_stochmap_TEST(Node & node,bool time);


	//get the scale
	double get_scale();
	//this should be in the form of the get_scale so (log(scale)*totalscale);
	void set_scale(double sc, double tsc);
	void set_run_with_scale(bool );

/*
	for timing things
 */
	double ti;
	double ti2;
	double ti3;
};

#endif /* BIOGEOTREE_H_ */
