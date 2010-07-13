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
#include <functional>
#include <numeric>

using namespace std;

#include "BioGeoTree.h"
#include "BioGeoTreeTools.h"
#include "BranchSegment.h"
#include "RateMatrixUtils.h"
#include "RateModel.h"
#include "AncSplit.h"

#include <Phyl/TreeTemplate.h>
#include <Phyl/TreeTemplateTools.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Utils/BppVector.h>
using namespace bpp;

namespace {
	inline double MAX(const double &a, const double &b)
	        {return b > a ? (b) : double(a);}
}

BioGeoTree::BioGeoTree(TreeTemplate<Node> * tr, vector<double> ps){
	seg = "segments";
	age = "age";
	dc = "dist_conditionals";
	en = "excluded_dists";
	andc = "anc_dist_conditionals";
	/*
        reverse bit
	 */
	revB  = "revB";
	/*
	end of the reverse bits
	*/
	store_p_matrices = false;
	use_stored_matrices = false;
	tree = tr;
	numofnodes = tree->getNumberOfNodes();
	periods = ps;
	TreeTemplateTools * ttt = new TreeTemplateTools();
	/*
	 * initialize each node with segments
	 */
	cout << "initializing nodes..." << endl;
	for (int i=0;i<numofnodes;i++){
		tree_get_node_from_id[i] = tree->getNode(i);
		Vector<BranchSegment> * segs = new Vector<BranchSegment>();
		tree_get_node_from_id[i]->setNodeProperty(seg,*segs);
		//tree->setNodeProperty(i,seg,*segs);
		Vector<vector<int> > * ens = new Vector<vector<int> >();
		//tree->setNodeProperty(i,en,*ens);
		tree_get_node_from_id[i]->setNodeProperty(en,*ens);
	}

	/*
	 * initialize the actual branch segments for each node
	 */
	cout << "initializing branch segments..." << endl;
	for (int i=0;i<numofnodes;i++){
		//if (tree->getNode(i)->hasFather()){
		if (tree_get_node_from_id[i]->hasFather()){
			vector<double> pers(periods);
			//double anc = ttt->getHeight(*tree->getNode(i)->getFather());
			//double des = ttt->getHeight(*tree->getNode(i));
			double anc = ttt->getHeight(*tree_get_node_from_id[i]->getFather());
			double des = ttt->getHeight(*tree_get_node_from_id[i]);
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
							//((bpp::Vector<BranchSegment>*) tree->getNodeProperty(i,seg))->push_back(tseg);
							((bpp::Vector<BranchSegment>*) tree_get_node_from_id[i]->getNodeProperty(seg))->push_back(tseg);
						}
						//t += pers[j];
						t += duration; // TODO: make sure that this is all working
					}
					if (t > anc || pers[j] > t){
						break;
					}
				}
			}else{
				//BranchSegment tseg = BranchSegment(tree->getNode(i)->getDistanceToFather(),0);
				//((bpp::Vector<BranchSegment>*) tree->getNodeProperty(i,seg))->push_back(tseg);
				BranchSegment tseg = BranchSegment(tree_get_node_from_id[i]->getDistanceToFather(),0);
				((bpp::Vector<BranchSegment>*) tree_get_node_from_id[i]->getNodeProperty(seg))->push_back(tseg);
			}
		}
	}
	delete ttt;
}

void BioGeoTree::set_store_p_matrices(bool i){
	store_p_matrices = i;
}

void BioGeoTree::set_use_stored_matrices(bool i){
	use_stored_matrices = i;
}

void BioGeoTree::set_default_model(RateModel * mod){
	rootratemodel = mod;
	for(int i=0;i<numofnodes;i++){
		//bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(i)->getNodeProperty(seg));
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree_get_node_from_id[i]->getNodeProperty(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
			Vector<double> * distconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).distconds = distconds;
			Vector<double> * ancdistconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
			tsegs->at(j).ancdistconds = ancdistconds;
        }
	}
	Vector<double> * distconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
	tree->getRootNode()->setNodeProperty(dc,*distconds);
	Vector<double> * ancdistconds = new Vector<double> (rootratemodel->getDists()->size(), 0);
	tree->getRootNode()->setNodeProperty(andc,*ancdistconds);
}

void BioGeoTree::update_default_model(RateModel * mod){
	rootratemodel = mod;
	for(int i=0;i<numofnodes;i++){
		//bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(i)->getNodeProperty(seg));
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree_get_node_from_id[i]->getNodeProperty(seg));
		for(unsigned int j=0;j<tsegs->size();j++){
			tsegs->at(j).setModel(mod);
		}
	}
}

void BioGeoTree::set_tip_conditionals(map<string,vector<int> > distrib_data){
	int numofleaves = tree->getNumberOfLeaves();
	vector<Node *> lvs = tree->getLeaves();
	for(int i=0;i<numofleaves;i++){
		//bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getLeaves().at(i)->getNodeProperty(seg));
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) lvs.at(i)->getNodeProperty(seg));
		RateModel * mod = tsegs->at(0).getModel();
		//int ind1 = get_vector_int_index_from_multi_vector_int(
		//		&distrib_data[tree->getLeaves().at(i)->getName()],mod->getDists());
		int ind1 = get_vector_int_index_from_multi_vector_int(
						&distrib_data[lvs.at(i)->getName()],mod->getDists());
		tsegs->at(0).distconds->at(ind1) = 1.0;
	}
}

void BioGeoTree::set_excluded_dist(vector<int> ind,Node * node){
	((bpp::Vector<vector<int> >*) node->getNodeProperty(en))->push_back(ind);
}

/*
 * **************************************************
 *
 *
 *
 * **************************************************
 */

double BioGeoTree::eval_likelihood(bool marginal){
	if( rootratemodel->sparse == true){
		columns = new vector<int>(rootratemodel->getDists()->size());
		whichcolumns = new vector<int>();
	}
	ancdist_conditional_lh(*tree->getRootNode(),marginal);
	if( rootratemodel->sparse == true){
		delete columns;
		delete whichcolumns;
	}
	return calculate_vector_double_sum(*
			(bpp::Vector<double>*) tree->getRootNode()->getNodeProperty(dc));

}

Vector<double> BioGeoTree::conditionals(Node & node, bool marginal,bool sparse){
	Vector<double> distconds;
	bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
	distconds = *tsegs->at(0).distconds;
	for(unsigned int i=0;i<tsegs->size();i++){
		for(unsigned int j=0;j<distconds.size();j++){
				tsegs->at(i).distconds->at(j) = distconds.at(j);
		}
		RateModel * rm = tsegs->at(i).getModel();
		Vector<double> * v = new Vector<double> (rootratemodel->getDists()->size(), 0);
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
					double maxnum = 0;
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

void BioGeoTree::ancdist_conditional_lh(Node & node, bool marginal){
	Vector<double> distconds(rootratemodel->getDists()->size(), 0);
	if (node.isLeaf()==false){//is not a tip
		Node * c1 = node.getSon(0);
		Node * c2 = node.getSon(1);
		RateModel * model;
		if(node.hasFather()==true){
			bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
			model = tsegs->at(0).getModel();
		}else{
			model = rootratemodel;
		}
		ancdist_conditional_lh(*c1,marginal);
		ancdist_conditional_lh(*c2,marginal);
		bool sparse = rootratemodel->sparse;
		Vector<double> v1;
		Vector<double> v2;
		if(sparse == true){
			//getcolumns
			bpp::Vector<BranchSegment>* c1tsegs = ((bpp::Vector<BranchSegment>*) c1->getNodeProperty(seg));
			bpp::Vector<BranchSegment>* c2tsegs = ((bpp::Vector<BranchSegment>*) c2->getNodeProperty(seg));
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
		v1 =conditionals(*c1,marginal,sparse);
		v2 =conditionals(*c2,marginal,sparse);

		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		//cl1 = clock();
		for (unsigned int i=0;i<dists->size();i++){
			if(accumulate(dists->at(i).begin(),dists->at(i).end(),0) > 0){
				double lh = 0.0;
				bpp::Vector<vector<int> >* exdist = ((bpp::Vector<vector<int> >*) node.getNodeProperty(en));
				int cou = count(exdist->begin(),exdist->end(),dists->at(i));
				if(cou == 0){
					iter_ancsplits_just_int(rootratemodel,dists->at(i),leftdists,rightdists,weight);
					for (unsigned int j=0;j<leftdists.size();j++){
						int ind1 = leftdists[j];
						int ind2 = rightdists[j];
						double lh_part = v1.at(ind1)*v2.at(ind2);
						lh += (lh_part * weight);
					}
				}
				distconds.at(i)= lh;
			}
		}
		///cl2 = clock();
		//ti += cl2-cl1;
	}else{
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
		distconds = *tsegs->at(0).distconds;
	}
	if(node.hasFather() == true){
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
		for(unsigned int i=0;i<distconds.size();i++){
			tsegs->at(0).distconds->at(i) = distconds.at(i);
		}
	}else{
		for(unsigned int i=0;i<distconds.size();i++){
			((Vector<double>*)node.getNodeProperty(dc))->at(i) = distconds.at(i);
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
void BioGeoTree::setFossilatNodeByMRCA(vector<string> nodeNames, int fossilarea){
	BioGeoTreeTools tt;
	vector<int> nodeIds;
	for(unsigned int i=0;i<nodeNames.size();i++){
		nodeIds.push_back(tree->getNode(nodeNames[i])->getId());
	}
	int id = tt.getLastCommonAncestor(*tree,nodeIds);
	vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			bpp::Vector<vector<int> > * exd = ((bpp::Vector<vector<int> > *) tree->getNodeProperty(id,en));
			exd->push_back(dists->at(i));
		}
	}
}
void BioGeoTree::setFossilatNodeByMRCA_id(int id, int fossilarea){
	vector<vector<int> > * dists = rootratemodel->getDists();
	for(unsigned int i=0;i<dists->size();i++){
		if(dists->at(i).at(fossilarea) == 0){
			bpp::Vector<vector<int> > * exd = ((bpp::Vector<vector<int> > *) tree->getNodeProperty(id,en));
			exd->push_back(dists->at(i));
		}
	}
}
void BioGeoTree::setFossilatBranchByMRCA(vector<string> nodeNames, int fossilarea, double age){
	BioGeoTreeTools tt;
	TreeTemplateTools * ttt = new TreeTemplateTools();
	vector<int> nodeIds;
	for(unsigned int i=0;i<nodeNames.size();i++){
		nodeIds.push_back(tree->getNode(nodeNames[i])->getId());
	}
	int id = tt.getLastCommonAncestor(*tree,nodeIds);
	//bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(id)->getNodeProperty(seg));
	//double startage = ttt->getHeight(*tree->getNode(id));
	bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree_get_node_from_id[id]->getNodeProperty(seg));
	double startage = ttt->getHeight(*tree_get_node_from_id[id]);
	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
	delete ttt;
}
void BioGeoTree::setFossilatBranchByMRCA_id(int id, int fossilarea, double age){
	TreeTemplateTools * ttt = new TreeTemplateTools();
	//bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree->getNode(id)->getNodeProperty(seg));
	//double startage = ttt->getHeight(*tree->getNode(id));
	bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) tree_get_node_from_id[id]->getNodeProperty(seg));
	double startage = ttt->getHeight(*tree_get_node_from_id[id]);

	for(unsigned int i=0;i<tsegs->size();i++){
		if(age > startage && age < (startage+tsegs->at(i).getDuration())){
			tsegs->at(i).setFossilArea(fossilarea);
		}
		startage += tsegs->at(i).getDuration();
	}
	delete ttt;
}


/************************************************************
 forward and reverse stuff
 ************************************************************/
//add joint
void BioGeoTree::prepare_ancstate_reverse(){
    reverse(*tree->getRootNode());
}

/*
 * called from prepare_ancstate_reverse and that is all
 */
void BioGeoTree::reverse(Node & node){
	Vector<double> * revconds = new Vector<double> (rootratemodel->getDists()->size(), 0);//need to delete this at some point
	if (node == *tree->getRootNode()) {
		for(unsigned int i=0;i<rootratemodel->getDists()->size();i++){
			revconds->at(i) = 1.0;//prior
		}
		node.setNodeProperty(revB,*revconds);
		for(unsigned int i = 0;i<node.getNumberOfSons();i++){
			reverse(*node.getSon(i));
		}
	}else if(node.isLeaf() == false){
		//calculate A i 
		//sum over all alpha k of sister node of the parent times the priors of the speciations 
		//(weights) times B of parent j
		Vector<double> * parrev = ((Vector<double>*)node.getFather()->getNodeProperty(revB));
		Vector<double> sisdistconds;
		if(node.getFather()->getSon(0) != &node){
			bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getFather()->getSon(0)->getNodeProperty(seg));
			sisdistconds = tsegs->at(0).alphas;
		}else{
			bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getFather()->getSon(1)->getNodeProperty(seg));
			sisdistconds = tsegs->at(0).alphas;
		}
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		//cl1 = clock();
		Vector<double> tempA (rootratemodel->getDists()->size(),0);
		for (unsigned int i = 0; i < dists->size(); i++) {
			if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
				bpp::Vector<vector<int> >* exdist =
						((bpp::Vector<vector<int> >*) node.getNodeProperty(en));
				int cou = count(exdist->begin(), exdist->end(), dists->at(i));
				if (cou == 0) {
					iter_ancsplits_just_int(rootratemodel, dists->at(i),
							leftdists, rightdists, weight);
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
		bpp::Vector<BranchSegment>* tsegs = ((bpp::Vector<BranchSegment>*) node.getNodeProperty(seg));
		for(unsigned int k=0;k<tsegs->size();k++){
			RateModel * rm = tsegs->at(k).getModel();
			vector<vector<double > > * p = &rm->stored_p_matrices[tsegs->at(k).getPeriod()][tsegs->at(k).getDuration()];
			for(unsigned int j=0;j < dists->size();j++){
				if(accumulate(dists->at(j).begin(), dists->at(j).end(), 0) > 0){
					for (unsigned int i = 0; i < dists->size(); i++) {
						if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
							revconds->at(j) += tempA[i]*((*p)[i][j]);
						}
					}
				}
			}
		}
		node.setNodeProperty(revB,*revconds);
		for(unsigned int i = 0;i<node.getNumberOfSons();i++){
			reverse(*node.getSon(i));
		}
	}//there should be no else
}

/*
 * calculates the most likely split (not state) -- the traditional result for lagrange
 */

map<vector<int>,vector<AncSplit> > BioGeoTree::calculate_ancsplit_reverse(Node & node,bool marg){
	Vector<double> * Bs = (Vector<double> *) node.getNodeProperty(revB);
	map<vector<int>,vector<AncSplit> > ret;
	for(unsigned int j=0;j<rootratemodel->getDists()->size();j++){
		vector<int> dist = rootratemodel->getDists()->at(j);
		vector<AncSplit> ans = iter_ancsplits(rootratemodel,dist);
		if (node.isLeaf()==false){//is not a tip
			Node * c1 = node.getSon(0);
			Node * c2 = node.getSon(1);
			bpp::Vector<BranchSegment>* tsegs1 = ((bpp::Vector<BranchSegment>*) c1->getNodeProperty(seg));
			bpp::Vector<BranchSegment>* tsegs2 = ((bpp::Vector<BranchSegment>*) c2->getNodeProperty(seg));
			for (unsigned int i=0;i<ans.size();i++){
				Vector<double> v1  =tsegs1->at(0).alphas;
				Vector<double> v2 = tsegs2->at(0).alphas;
				double lh = (v1[ans[i].ldescdistint]*v2[ans[i].rdescdistint]*Bs->at(j)*ans[i].getWeight());
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

vector<double> BioGeoTree::calculate_ancstate_reverse(Node & node,bool marg){
	if (node.isLeaf()==false){//is not a tip
		Vector<double> * Bs = (Vector<double> *) node.getNodeProperty(revB);
		vector<vector<int> > * dists = rootratemodel->getDists();
		vector<int> leftdists;
		vector<int> rightdists;
		double weight;
		Node * c1 = node.getSon(0);
		Node * c2 = node.getSon(1);
		bpp::Vector<BranchSegment>* tsegs1 = ((bpp::Vector<BranchSegment>*) c1->getNodeProperty(seg));
		bpp::Vector<BranchSegment>* tsegs2 = ((bpp::Vector<BranchSegment>*) c2->getNodeProperty(seg));
		Vector<double> v1  =tsegs1->at(0).alphas;
		Vector<double> v2 = tsegs2->at(0).alphas;
		//cl1 = clock();
		Vector<double> LHOODS (rootratemodel->getDists()->size(),0);
		for (unsigned int i = 0; i < dists->size(); i++) {
			if (accumulate(dists->at(i).begin(), dists->at(i).end(), 0) > 0) {
				bpp::Vector<vector<int> >* exdist =
				((bpp::Vector<vector<int> >*) node.getNodeProperty(en));
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


