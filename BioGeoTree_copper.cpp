/*
 * BioGeoTree.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */
#include <vector>
#include <string>
#include <algorithm>
#include <ctime>
//#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
using namespace std;

#include <armadillo>
using namespace arma;

#include "BioGeoTree_copper.h"
#include "BioGeoTreeTools_copper.h"
#include "BranchSegment_copper.h"
#include "RateMatrixUtils.h"
#include "RateModel.h"
#include "AncSplit.h"

#include "tree.h"
#include "node.h"
#include "vector_node_object.h"

//octave usage
//#include <octave/oct.h>

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif


#ifdef BIGTREE
namespace{
	inline mpfr_class MAX(const mpfr_class &a, const mpfr_class &b)
	        {return b > a ? (b) : mpfr_class(a);}
}
#else
namespace {
	inline double MAX(const double &a, const double &b)
	        {return b > a ? (b) : double(a);}
}
#endif

/*
 * sloppy beginning but best for now because of the complicated bits
 */

BioGeoTree_copper::BioGeoTree_copper(Tree * tr, vector<double> ps):tree(tr),periods(ps),
		seg("segments"),age("age"),dc("dist_conditionals"),en("excluded_dists"),
		andc("anc_dist_conditionals"),columns(NULL),whichcolumns(NULL),rootratemodel(NULL),
		distmap(NULL),store_p_matrices(false),use_stored_matrices(false),revB("revB"),
		rev(false),rev_exp_number("rev_exp_number"),rev_exp_time("rev_exp_time"),
		stochastic(false),stored_EN_matrices(map<int,map<double, mat > >()),
		stored_ER_matrices(map<int,map<double, mat > >()){

	/*
	 * initialize each node with segments
	 */
	cout << "initializing nodes..." << endl;
	for(int i=0;i<tree->getNodeCount();i++){
		if(tree->getNode(i)->getBL()<0.000001)
			tree->getNode(i)->setBL(0.000001);
		VectorNodeObject<BranchSegment> * segs = new VectorNodeObject<BranchSegment>();
		tree->getNode(i)->assocObject(seg,*segs);
		delete segs;
		VectorNodeObject<vector<int> > * ens = new VectorNodeObject<vector<int> >();
		tree->getNode(i)->assocObject(en,*ens);
		delete ens;
	}
	/*
	 * initialize the actual branch segments for each node
	 */
	tree->setHeightFromTipToNodes();
	cout << "initializing branch segments..." << endl;
	for(int i=0;i<tree->getNodeCount();i++){
		if (tree->getNode(i)->hasParent()){
			vector<double> pers(periods);
			double anc = tree->getNode(i)->getParent()->getHeight();
			double des = tree->getNode(i)->getHeight();
			//assert anc > des:q
			double t = des;
			if (pers.size() > 0){
				for(unsigned int j=0;j<pers.size();j++){
					double s = 0;
					if(pers.size() == 1)
						s = pers[0];
					for (unsigned int k=0;k<j+1;k++){
						s += pers[k];
					}
					if (t < s){
						double duration = min(s-t,anc-t);
						if (duration > 0){
							BranchSegment tseg = BranchSegment(duration,j);
							((VectorNodeObject<BranchSegment>*) tree->getNode(i)->getObject(seg))->push_back(tseg);
						}
						t += duration; // TODO: make sure that this is all working
					}
					if (t > anc || pers[j] > t){
						break;
					}
				}
			}else{
				BranchSegment tseg = BranchSegment(tree->getNode(i)->getBL(),0);
				((VectorNodeObject<BranchSegment>*) tree->getNode(i)->getObject(seg))->push_back(tseg);
			}
		}
	}
}

void BioGeoTree_copper::set_store_p_matrices(bool i){
	store_p_matrices = i;
}

void BioGeoTree_copper::set_use_stored_matrices(bool i){
	use_stored_matrices = i;
}

void BioGeoTree_copper::set_default_model(RateModel * mod){
	rootratemodel = mod;
	for(int i=0;i<tree->getNodeCount();i++){
			VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getNode(i)->getObject(seg));
			for(unsigned int j=0;j<tsegs->size();j++){
				tsegs->at(j).setModel(mod);
#ifdef BIGTREE
				VectorNodeObject<mpfr_class> * distconds = new VectorNodeObject<mpfr_class> (rootratemodel->getDists()->size(), 0);
				tsegs->at(j).distconds = distconds;
				VectorNodeObject<mpfr_class> * ancdistconds = new VectorNodeObject<mpfr_class> (rootratemodel->getDists()->size(), 0);
				tsegs->at(j).ancdistconds = ancdistconds;
#else
				VectorNodeObject<double> * distconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
				tsegs->at(j).distconds = distconds;
				VectorNodeObject<double> * ancdistconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
				tsegs->at(j).ancdistconds = ancdistconds;
#endif
			}
		}
#ifdef BIGTREE
	VectorNodeObject<mpfr_class> * distconds = new VectorNodeObject<mpfr_class> (rootratemodel->getDists()->size(), 0);
	tree->getRoot()->assocObject(dc,*distconds);
	delete distconds;
	VectorNodeObject<mpfr_class> * ancdistconds = new VectorNodeObject<mpfr_class> (rootratemodel->getDists()->size(), 0);
	tree->getRoot()->assocObject(andc,*ancdistconds);
	delete ancdistconds;
#else
	VectorNodeObject<double> * distconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
	tree->getRoot()->assocObject(dc,*distconds);
	delete distconds;
	VectorNodeObject<double> * ancdistconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
	tree->getRoot()->assocObject(andc,*ancdistconds);
	delete ancdistconds;
#endif
}

void BioGeoTree_copper::update_default_model(RateModel * mod){
	rootratemodel = mod;

	for(int i=0;i<tree->getNodeCount();i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getNode(i)->getObject(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
		}
	}
}

void BioGeoTree_copper::set_tip_conditionals(map<string,vector<int> > distrib_data){
	int numofleaves = tree->getExternalNodeCount();
	for(int i=0;i<numofleaves;i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getExternalNode(i)->getObject(seg));
		RateModel * mod = tsegs->at(0).getModel();
		int ind1 = get_vector_int_index_from_multi_vector_int(
						&distrib_data[tree->getExternalNode(i)->getName()],mod->getDists());
		tsegs->at(0).distconds->at(ind1) = 1.0;
	}
}

void BioGeoTree_copper::set_excluded_dist(vector<int> ind,Node * node){
	((VectorNodeObject<vector<int> >*) node->getObject(en))->push_back(ind);
}

/*
 * **************************************************
 *
 *
 *
 * **************************************************
 */
#ifdef BIGTREE
mpfr_class calculate_vector_mpfr_class_sum(vector<mpfr_class> & in);
mpfr_class calculate_vector_mpfr_class_sum(vector<mpfr_class> & in){
	mpfr_class sum = 0;
	for (unsigned int i=0;i<in.size();i++){
		sum += in[i];
	}
	return sum;
}
#endif

double BioGeoTree_copper::eval_likelihood(bool marginal){
	if( rootratemodel->sparse == true){
		columns = new vector<int>(rootratemodel->getDists()->size());
		whichcolumns = new vector<int>();
	}
	ancdist_conditional_lh(*tree->getRoot(),marginal);
	if( rootratemodel->sparse == true){
		delete columns;
		delete whichcolumns;
	}
#ifdef BIGTREE
	mpfr_class f = 	(-log(calculate_vector_mpfr_class_sum(*
			(VectorNodeObject<mpfr_class>*) tree->getRoot()->getObject(dc))));
	double x = f.get_d();
	return x;
#else
	return (-log(calculate_vector_double_sum(*
			(VectorNodeObject<double>*) tree->getRoot()->getObject(dc))));
#endif
}

#ifdef BIGTREE
VectorNodeObject<mpfr_class> BioGeoTree_copper::conditionals(Node & node, bool marginal,bool sparse){
#else
VectorNodeObject<double> BioGeoTree_copper::conditionals(Node & node, bool marginal,bool sparse){
#endif
#ifdef BIGTREE
	VectorNodeObject<mpfr_class> distconds;
#else
	VectorNodeObject<double> distconds;
#endif
	VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
	distconds = *tsegs->at(0).distconds;
	for(unsigned int i=0;i<tsegs->size();i++){
		for(unsigned int j=0;j<distconds.size();j++){
				tsegs->at(i).distconds->at(j) = distconds.at(j);
		}
		RateModel * rm = tsegs->at(i).getModel();
#ifdef BIGTREE
		VectorNodeObject<mpfr_class> * v = new VectorNodeObject<mpfr_class> (rootratemodel->getDists()->size(), 0);
#else
		VectorNodeObject<double> * v = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);
#endif
		vector<int> distrange;
		if(tsegs->at(i).get_start_dist_int() != -666){
			int ind1 = tsegs->at(i).get_start_dist_int();
			distrange.push_back(ind1);
		}else if(tsegs->at(i).getFossilAreas().size()>0){
			for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
				distrange.push_back(j);
			}
			for(unsigned int k=0;k<distrange.size();k++){
				bool flag = true;
				for(unsigned int x = 0;x<tsegs->at(i).getFossilAreas().size();x++){
					if (tsegs->at(i).getFossilAreas()[x] == 1 && distrange.at(x) == 0){
						flag = false;
					}
				}
				if(flag == true){
                    distrange.erase(distrange.begin()+k);
				}
			}
		}else{
			for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
				distrange.push_back(j);
			}
		}
		/*
		 * marginal
		 */
		if(marginal == true){
			if(sparse == false){
				vector<vector<double > > p;
				if(use_stored_matrices == false){
					p= rm->setup_fortran_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),
																 store_p_matrices);
				}else{
					p = rm->stored_p_matrices[tsegs->at(i).getPeriod()][tsegs->at(i).getDuration()];
				}
				for(unsigned int j=0;j<distrange.size();j++){
					for(unsigned int k=0;k<distconds.size();k++){
						v->at(distrange[j]) += (distconds.at(k)*p[distrange[j]][k]);
					}
				}
			}else{//sparse
				/*
				 testing pthread version
				 */
				if(rm->get_nthreads() > 0){
					vector<vector<double > > p = rm->setup_pthread_sparse_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),*whichcolumns);
					for(unsigned int j=0;j<distrange.size();j++){
						for(unsigned int k=0;k<distconds.size();k++){
							v->at(distrange[j]) += (distconds.at(k)*p[distrange[j]][k]);
						}
					}
				}else{
					for(unsigned int j=0;j<distrange.size();j++){
						bool inthere = false;
						if(columns->at(distrange[j]) == 1)
							inthere = true;
						vector<double > p;
						if(inthere == true){
							p = rm->setup_sparse_single_column_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),distrange[j]);
						}else{
							p = vector<double>(distconds.size(),0);
						}
						for(unsigned int k=0;k<distconds.size();k++){
							v->at(distrange[j]) += (distconds.at(k)*p[k]);
						}
					}
				}
			}
		}
		/*
		 * joint reconstruction
		 * NOT FINISHED YET -- DONT USE
		 */
		else{
			if(sparse == false){
				vector<vector<double > > p = rm->setup_fortran_P(tsegs->at(i).getPeriod(),tsegs->at(i).getDuration(),store_p_matrices);
				for(unsigned int j=0;j<distrange.size();j++){
#ifdef BIGTREE
					mpfr_class maxnum = 0;
#else
					double maxnum = 0;
#endif
					for(unsigned int k=0;k<distconds.size();k++){
						maxnum = MAX((distconds.at(k)*p[distrange[j]][k]),maxnum);
					}
					v->at(distrange[j]) = maxnum;
				}
			}else{//sparse

			}
		}
		for(unsigned int j=0;j<distconds.size();j++){
			distconds[j] = v->at(j);
		}
		if(store_p_matrices == true){
			tsegs->at(i).seg_sp_alphas = distconds;
		}
		delete v;
	}
	/*
	 * if store is true we want to store the conditionals for each node
	 * for possible use in ancestral state reconstruction
	 */
	if(store_p_matrices == true){
		tsegs->at(0).alphas = distconds;
	}
	return distconds;
}

void BioGeoTree_copper::ancdist_conditional_lh(Node & node, bool marginal){
#ifdef BIGTREE
	VectorNodeObject<mpfr_class> distconds(rootratemodel->getDists()->size(), 0);
#else
	VectorNodeObject<double> distconds(rootratemodel->getDists()->size(), 0);
#endif
	if (node.isExternal()==false){//is not a tip
		Node * c1 = &node.getChild(0);
		Node * c2 = &node.getChild(1);
		RateModel * model;
		if(node.hasParent()==true){
			VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
			model = tsegs->at(0).getModel();
		}else{
			model = rootratemodel;
		}
		ancdist_conditional_lh(*c1,marginal);
		ancdist_conditional_lh(*c2,marginal);
		bool sparse = rootratemodel->sparse;
#ifdef BIGTREE
		VectorNodeObject<mpfr_class> v1;
		VectorNodeObject<mpfr_class> v2;
#else
		VectorNodeObject<double> v1;
		VectorNodeObject<double> v2;
#endif
#ifdef BIGTREE

#else
		if(sparse == true){
			//getcolumns
			VectorNodeObject<BranchSegment>* c1tsegs = ((VectorNodeObject<BranchSegment>*) c1->getObject(seg));
			VectorNodeObject<BranchSegment>* c2tsegs = ((VectorNodeObject<BranchSegment>*) c2->getObject(seg));
			vector<int> lcols = get_columns_for_sparse(*c1tsegs->at(0).distconds,rootratemodel);
			vector<int> rcols = get_columns_for_sparse(*c2tsegs->at(0).distconds,rootratemodel);
			whichcolumns->clear();
			for(unsigned int i=0;i<lcols.size();i++){
				if(lcols[i]==1 || rcols[i] ==1){
					columns->at(i)=1;
					if(i!=0 && count(whichcolumns->begin(),whichcolumns->end(),i) == 0)
						whichcolumns->push_back(i);
				}else{
					columns->at(i)=0;
				}
			}
			if(calculate_vector_int_sum(columns)==0){
				for(unsigned int i=0;i<lcols.size();i++){
					columns->at(i)=1;
				}
			}
			columns->at(0) = 0;
		}
#endif
		v1 =conditionals(*c1,marginal,sparse);
		v2 =conditionals(*c2,marginal,sparse);

		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		//cl1 = clock();
		for (unsigned int i=0;i<dists->size();i++){
			if(accumulate(dists->at(i).begin(),dists->at(i).end(),0) > 0){
#ifdef BIGTREE
				mpfr_class lh = 0.0;
#else
				double lh = 0.0;
#endif
				VectorNodeObject<vector<int> >* exdist = ((VectorNodeObject<vector<int> >*) node.getObject(en));
				int cou = count(exdist->begin(),exdist->end(),dists->at(i));
				if(cou == 0){
					iter_ancsplits_just_int(rootratemodel,dists->at(i),leftdists,rightdists,weight);
					for (unsigned int j=0;j<leftdists.size();j++){
						int ind1 = leftdists[j];
						int ind2 = rightdists[j];
#ifdef BIGTREE
						mpfr_class lh_part = v1.at(ind1)*v2.at(ind2);
#else
						double lh_part = v1.at(ind1)*v2.at(ind2);
#endif
						lh += (lh_part * weight);
					}
				}
				distconds.at(i)= lh;
			}
		}
		///cl2 = clock();
		//ti += cl2-cl1;
	}else{
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		distconds = *tsegs->at(0).distconds;
	}
	if(node.hasParent() == true){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		for(unsigned int i=0;i<distconds.size();i++){
			tsegs->at(0).distconds->at(i) = distconds.at(i);
		}
	}else{
		for(unsigned int i=0;i<distconds.size();i++){
#ifdef BIGTREE
			((VectorNodeObject<mpfr_class>*)node.getObject(dc))->at(i) = distconds.at(i);
#else
			((VectorNodeObject<double>*)node.getObject(dc))->at(i) = distconds.at(i);
#endif
		}
	}
}

/*
 * ********************************************
 *
 * adds fossils either at the node or along a branch
 *
 * ********************************************
 */
void BioGeoTree_copper::setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea){
	Node * mrca = tree->getMRCA(nodeNames);
	vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			VectorNodeObject<vector<int> > * exd = ((VectorNodeObject<vector<int> > *) mrca->getObject(en));
			exd->push_back(dists->at(i));
		}
	}
}
void BioGeoTree_copper::setFossilatNodeByMRCA_id(Node * id, int fossilarea){
	vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			VectorNodeObject<vector<int> > * exd = ((VectorNodeObject<vector<int> > *) id->getObject(en));
			exd->push_back(dists->at(i));
		}
	}
}
void BioGeoTree_copper::setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea, double age){
	Node * mrca = tree->getMRCA(nodeNames);
	VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) mrca->getObject(seg));
	double startage = mrca->getHeight();
	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
}
void BioGeoTree_copper::setFossilatBranchByMRCA_id(Node * id, int fossilarea, double age){
	VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) id->getObject(seg));
	double startage = id->getHeight();

	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
}


/************************************************************
 forward and reverse stuff for ancestral states
 ************************************************************/
//add joint
void BioGeoTree_copper::prepare_ancstate_reverse(){
    reverse(*tree->getRoot());
}

/*
 * called from prepare_ancstate_reverse and that is all
 */
void BioGeoTree_copper::reverse(Node & node){
	rev = true;
#ifdef BIGTREE
	VectorNodeObject<mpfr_class> * revconds = new VectorNodeObject<mpfr_class> (rootratemodel->getDists()->size(), 0);//need to delete this at some point
#else
	VectorNodeObject<double> * revconds = new VectorNodeObject<double> (rootratemodel->getDists()->size(), 0);//need to delete this at some point
#endif
	if (&node == tree->getRoot()) {
		for(unsigned int i=0;i<rootratemodel->getDists()->size();i++){
			revconds->at(i) = 1.0;//prior
		}
		node.assocObject(revB,*revconds);
		delete revconds;
		for(int i = 0;i<node.getChildCount();i++){
			reverse(node.getChild(i));
		}
	}else{
	//else if(node.isExternal() == false){
		//calculate A i 
		//sum over all alpha k of sister node of the parent times the priors of the speciations 
		//(weights) times B of parent j
#ifdef BIGTREE
		VectorNodeObject<mpfr_class> * parrev = ((VectorNodeObject<mpfr_class>*)node.getParent()->getObject(revB));
		VectorNodeObject<mpfr_class> sisdistconds;
#else
		VectorNodeObject<double> * parrev = ((VectorNodeObject<double>*)node.getParent()->getObject(revB));
		VectorNodeObject<double> sisdistconds;
#endif
		if(&node.getParent()->getChild(0) != &node){
			VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getParent()->getChild(0).getObject(seg));
			sisdistconds = tsegs->at(0).alphas;
		}else{
			VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getParent()->getChild(1).getObject(seg));
			sisdistconds = tsegs->at(0).alphas;
		}
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		//cl1 = clock();
#ifdef BIGTREE
		VectorNodeObject<mpfr_class> tempA (rootratemodel->getDists()->size(),0);
#else
		VectorNodeObject<double> tempA (rootratemodel->getDists()->size(),0);
#endif
		for (unsigned int i = 0; i < dists->size(); i++) {
			if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
				VectorNodeObject<vector<int> >* exdist = ((VectorNodeObject<vector<int> >*) node.getObject(en));
				int cou = count(exdist->begin(), exdist->end(), dists->at(i));
				if (cou == 0) {
					iter_ancsplits_just_int(rootratemodel, dists->at(i), leftdists, rightdists, weight);
					//root has i, curnode has left, sister of cur has right
					for (unsigned int j = 0; j < leftdists.size(); j++) {
						int ind1 = leftdists[j];
						int ind2 = rightdists[j];
						tempA[ind1] += (sisdistconds.at(ind2)*weight*parrev->at(i));
					}
				}
			}
		}

		//now calculate node B
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		vector<double> tempmoveA(tempA);
		//for(unsigned int ts=0;ts<tsegs->size();ts++){
		for(int ts = tsegs->size()-1;ts != -1;ts--){
			for(unsigned int j=0;j<dists->size();j++){revconds->at(j) = 0;}
			RateModel * rm = tsegs->at(ts).getModel();
			vector<vector<double > > * p = &rm->stored_p_matrices[tsegs->at(ts).getPeriod()][tsegs->at(ts).getDuration()];
			mat * EN = NULL;
			mat * ER = NULL;
			VectorNodeObject<double> tempmoveAer(tempA);
			VectorNodeObject<double> tempmoveAen(tempA);
			if(stochastic == true){
				//initialize the segment B's
				for(unsigned int j=0;j<dists->size();j++){tempmoveAer[j] = 0;}
				for(unsigned int j=0;j<dists->size();j++){tempmoveAen[j] = 0;}
				EN = &stored_EN_matrices[tsegs->at(ts).getPeriod()][tsegs->at(ts).getDuration()];
				ER = &stored_ER_matrices[tsegs->at(ts).getPeriod()][tsegs->at(ts).getDuration()];
			}
			for(unsigned int j=0;j < dists->size();j++){
				if(accumulate(dists->at(j).begin(), dists->at(j).end(), 0) > 0){
					for (unsigned int i = 0; i < dists->size(); i++) {
						if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
							revconds->at(j) += tempmoveA[i]*((*p)[i][j]);//tempA needs to change each time
							if(stochastic == true){
								tempmoveAer[j] += tempmoveA[i]*(((*ER)(i,j)));
								tempmoveAen[j] += tempmoveA[i]*(((*EN)(i,j)));
							}
						}
					}
				}
			}
			for(unsigned int j=0;j<dists->size();j++){tempmoveA[j] = revconds->at(j);}
			if(stochastic == true){
				tsegs->at(ts).seg_sp_stoch_map_revB_time = tempmoveAer;
				tsegs->at(ts).seg_sp_stoch_map_revB_number = tempmoveAen;
			}
		}
		node.assocObject(revB,*revconds);
		delete revconds;
		for(int i = 0;i<node.getChildCount();i++){
			reverse(node.getChild(i));
		}
	}
}

/*
 * calculates the most likely split (not state) -- the traditional result for lagrange
 */

map<vector<int>,vector<AncSplit> > BioGeoTree_copper::calculate_ancsplit_reverse(Node & node,bool marg){
#ifdef BIGTREE
	VectorNodeObject<mpfr_class> * Bs = (VectorNodeObject<mpfr_class> *) node.getObject(revB);
#else
	VectorNodeObject<double> * Bs = (VectorNodeObject<double> *) node.getObject(revB);
#endif
	map<vector<int>,vector<AncSplit> > ret;
	for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
		vector<int> dist = rootratemodel->getDists()->at(j);
		vector<AncSplit> ans = iter_ancsplits(rootratemodel,dist);
		if (node.isExternal()==false){//is not a tip
			Node * c1 = &node.getChild(0);
			Node * c2 = &node.getChild(1);
			VectorNodeObject<BranchSegment>* tsegs1 = ((VectorNodeObject<BranchSegment>*) c1->getObject(seg));
			VectorNodeObject<BranchSegment>* tsegs2 = ((VectorNodeObject<BranchSegment>*) c2->getObject(seg));
			for (unsigned int i=0;i<ans.size();i++){
#ifdef BIGTREE
				VectorNodeObject<mpfr_class> v1  =tsegs1->at(0).alphas;
				VectorNodeObject<mpfr_class> v2 = tsegs2->at(0).alphas;
				mpfr_class lh = (v1[ans[i].ldescdistint]*v2[ans[i].rdescdistint]*Bs->at(j)*ans[i].getWeight());
#else
				VectorNodeObject<double> v1  =tsegs1->at(0).alphas;
				VectorNodeObject<double> v2 = tsegs2->at(0).alphas;
				double lh = (v1[ans[i].ldescdistint]*v2[ans[i].rdescdistint]*Bs->at(j)*ans[i].getWeight());
#endif
				ans[i].setLikelihood(lh);
				//cout << lh << endl;
			}
		}
		ret[dist] = ans;
	}
	return ret;
}

/*
 * calculates the ancestral area over all the possible splits
 */
#ifdef BIGTREE
vector<mpfr_class> BioGeoTree_copper::calculate_ancstate_reverse(Node & node,bool marg)
#else
vector<double> BioGeoTree_copper::calculate_ancstate_reverse(Node & node,bool marg)
#endif
	{
	if (node.isExternal()==false){//is not a tip
#ifdef BIGTREE
		VectorNodeObject<mpfr_class> * Bs = (VectorNodeObject<mpfr_class> *) node.getObject(revB);
#else
		VectorNodeObject<double> * Bs = (VectorNodeObject<double> *) node.getObject(revB);
#endif
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		Node * c1 = &node.getChild(0);
		Node * c2 = &node.getChild(1);
		VectorNodeObject<BranchSegment>* tsegs1 = ((VectorNodeObject<BranchSegment>*) c1->getObject(seg));
		VectorNodeObject<BranchSegment>* tsegs2 = ((VectorNodeObject<BranchSegment>*) c2->getObject(seg));
#ifdef BIGTREE
		VectorNodeObject<mpfr_class> v1  =tsegs1->at(0).alphas;
		VectorNodeObject<mpfr_class> v2 = tsegs2->at(0).alphas;
		VectorNodeObject<mpfr_class> LHOODS (dists->size(),0);
#else
		VectorNodeObject<double> v1  =tsegs1->at(0).alphas;
		VectorNodeObject<double> v2 = tsegs2->at(0).alphas;
		VectorNodeObject<double> LHOODS (dists->size(),0);
#endif
		for (unsigned int i = 0; i < dists->size(); i++) {
			if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
				VectorNodeObject<vector<int> >* exdist =
				((VectorNodeObject<vector<int> >*) node.getObject(en));
				int cou = count(exdist->begin(), exdist->end(), dists->at(i));
				if (cou == 0) {
					iter_ancsplits_just_int(rootratemodel, dists->at(i),
											leftdists, rightdists, weight);
					for (unsigned int j=0;j<leftdists.size();j++){
						int ind1 = leftdists[j];
						int ind2 = rightdists[j];
						LHOODS[i] += (v1.at(ind1)*v2.at(ind2)*weight);
					}
					LHOODS[i] *= Bs->at(i);
				}
			}
		}
		return LHOODS;
	}
}

/**********************************************************
 * forward and reverse stuff for stochastic mapping
 **********************************************************/

void BioGeoTree_copper::prepare_stochmap_reverse_all_nodes(int from , int to){
	stochastic = true;
	int ndists = rootratemodel->getDists()->size();

	//calculate and store local expectation matrix for each branch length
	for(int k = 0; k < tree->getNodeCount(); k++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getNode(k)->getObject(seg));
		for (unsigned int l = 0;l<tsegs->size();l++){
			int per = (*tsegs)[l].getPeriod();
			double dur =  (*tsegs)[l].getDuration();
			cx_mat eigvec(ndists,ndists);eigvec.fill(0);
			cx_mat eigval(ndists,ndists);eigval.fill(0);
			bool isImag = rootratemodel->get_eigenvec_eigenval_from_Q(&eigval, &eigvec,per);
			mat Ql(ndists,ndists);Ql.fill(0);Ql(from,to) = rootratemodel->get_Q()[per][from][to];
			mat W(ndists,ndists);W.fill(0);W(from,from) = 1;
			cx_mat summed(ndists,ndists);summed.fill(0);
			cx_mat summedR(ndists,ndists);summedR.fill(0);
			for(int i=0;i<ndists;i++){
				mat Ei(ndists,ndists);Ei.fill(0);Ei(i,i)=1;
				cx_mat Si(ndists,ndists);
				Si = eigvec * Ei * inv(eigvec);
				for(int j=0;j<ndists;j++){
					cx_double dij = (eigval(i,i)-eigval(j,j)) * dur;
					mat Ej(ndists,ndists);Ej.fill(0);Ej(j,j)=1;
					cx_mat Sj(ndists,ndists);
					Sj = eigvec * Ej * inv(eigvec);
					cx_double Iijt = 0;
					if (abs(dij) > 10){
						Iijt = (exp(eigval(i,i)*dur)-exp(eigval(j,j)*dur))/(eigval(i,i)-eigval(j,j));
					}else if(abs(dij) < 10e-20){
						Iijt = dur*exp(eigval(j,j)*dur)*(1.+dij/2.+pow(dij,2.)/6.+pow(dij,3.)/24.);
					}else{
						if(eigval(i,i) == eigval(j,j)){
							//WAS Iijt = dur*exp(eigval(j,j)*dur)*expm1(dij)/dij;
							if (isImag)
								Iijt = dur*exp(eigval(j,j)*dur)*(exp(dij)-1.)/dij;
							else
								Iijt = dur*exp(eigval(j,j)*dur)*(expm1(real(dij)))/dij;
						}else{
							//WAS Iijt = -dur*exp(eigval(i,i)*dur)*expm1(-dij)/dij;
							if (isImag)
								Iijt = -dur*exp(eigval(i,i)*dur)*(exp(-dij)-1.)/dij;
							else
								Iijt = -dur*exp(eigval(i,i)*dur)*(expm1(real(-dij)))/dij;
						}
					}
					summed += (Si  * Ql * Sj * Iijt);
					summedR += (Si * W * Sj * Iijt);
				}
			}
			stored_EN_matrices[per][dur] = (real(summed));
			stored_ER_matrices[per][dur] = (real(summedR));
		}
	}
}

/*
 * called directly after reverse_stochastic
 */

vector<double> BioGeoTree_copper::calculate_reverse_stochmap(Node & node,bool time){
	if (node.isExternal()==false){//is not a tip
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<double> totalExp (dists->size(),0);
		for(int t = 0;t<tsegs->size();t++){
			if (t == 0){
				vector<double> Bs;
				if(time)
					Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
				else
					Bs =  tsegs->at(t).seg_sp_stoch_map_revB_number;
				vector<int> leftdists;
				vector<int> rightdists;
				double weight;
				Node * c1 = &node.getChild(0);
				Node * c2 = &node.getChild(1);
				VectorNodeObject<BranchSegment>* tsegs1 = ((VectorNodeObject<BranchSegment>*) c1->getObject(seg));
				VectorNodeObject<BranchSegment>* tsegs2 = ((VectorNodeObject<BranchSegment>*) c2->getObject(seg));
				VectorNodeObject<double> v1  =tsegs1->at(0).alphas;
				VectorNodeObject<double> v2 = tsegs2->at(0).alphas;
				VectorNodeObject<double> LHOODS (dists->size(),0);
				for (unsigned int i = 0; i < dists->size(); i++) {
					if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
						VectorNodeObject<vector<int> >* exdist =
						((VectorNodeObject<vector<int> >*) node.getObject(en));
						int cou = count(exdist->begin(), exdist->end(), dists->at(i));
						if (cou == 0) {
							iter_ancsplits_just_int(rootratemodel, dists->at(i),
													leftdists, rightdists, weight);
							for (unsigned int j=0;j<leftdists.size();j++){
								int ind1 = leftdists[j];
								int ind2 = rightdists[j];
								LHOODS[i] += (v1.at(ind1)*v2.at(ind2)*weight);
							}
							LHOODS[i] *= Bs.at(i);
						}
					}
				}
				for(int i=0;i<dists->size();i++){
					totalExp[i] = LHOODS[i];
				}
			}else{
				vector<double> alphs = tsegs->at(t-1).seg_sp_alphas;
				vector<double> Bs;
				if(time)
					Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
				else
					Bs =  tsegs->at(t).seg_sp_stoch_map_revB_number;
				VectorNodeObject<double> LHOODS (dists->size(),0);
				for (unsigned int i = 0; i < dists->size(); i++) {
					if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
						VectorNodeObject<vector<int> >* exdist =
								((VectorNodeObject<vector<int> >*) node.getObject(en));
						int cou = count(exdist->begin(), exdist->end(), dists->at(i));
						if (cou == 0) {
							LHOODS[i] = Bs.at(i) * (alphs[i] );//do i do this or do i do from i to j
						}
					}
				}
				for(int i=0;i<dists->size();i++){
					totalExp[i] += LHOODS[i];
				}
			}
		}

		return totalExp;
	}else{
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) node.getObject(seg));
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<double> totalExp (dists->size(),0);
		for(int t = 0;t<tsegs->size();t++){
			if(t == 0){
				vector<double> Bs;
				if(time)
					Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
				else
					Bs =  tsegs->at(t).seg_sp_stoch_map_revB_number;
				VectorNodeObject<double> LHOODS (dists->size(),0);
				for (unsigned int i = 0; i < dists->size(); i++) {
					if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
						VectorNodeObject<vector<int> >* exdist =
								((VectorNodeObject<vector<int> >*) node.getObject(en));
						int cou = count(exdist->begin(), exdist->end(), dists->at(i));
						if (cou == 0) {
							LHOODS[i] = Bs.at(i) * (tsegs->at(0).distconds->at(i) );
						}
					}
				}
				for(int i=0;i<dists->size();i++){
					totalExp[i] = LHOODS[i];
				}
			}else{
				vector<double> alphs = tsegs->at(t-1).seg_sp_alphas;
				vector<double> Bs;
				if(time)
					Bs = tsegs->at(t).seg_sp_stoch_map_revB_time;
				else
					Bs =  tsegs->at(t).seg_sp_stoch_map_revB_number;
				VectorNodeObject<double> LHOODS (dists->size(),0);
				for (unsigned int i = 0; i < dists->size(); i++) {
					if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
						VectorNodeObject<vector<int> >* exdist =
								((VectorNodeObject<vector<int> >*) node.getObject(en));
						int cou = count(exdist->begin(), exdist->end(), dists->at(i));
						if (cou == 0) {
							LHOODS[i] = Bs.at(i) * (alphs[i]);
						}
					}
				}
				for(int i=0;i<dists->size();i++){
					totalExp[i] += LHOODS[i];
				}
			}
		}
		return totalExp;
	}
}

/**********************************************************
 * trash collection
 **********************************************************/
BioGeoTree_copper::~BioGeoTree_copper(){
	for(int i=0;i<tree->getNodeCount();i++){
		VectorNodeObject<BranchSegment>* tsegs = ((VectorNodeObject<BranchSegment>*) tree->getNode(i)->getObject(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			delete tsegs->at(j).distconds;
			delete tsegs->at(j).ancdistconds;
		}
		delete tree->getNode(i)->getObject(seg);
		delete tree->getNode(i)->getObject(en);
		if(rev == true && tree->getNode(i)->isInternal()){
			delete tree->getNode(i)->getObject(revB);
		}
	}
	delete tree->getRoot()->getObject(dc);
	delete tree->getRoot()->getObject(andc);
}

