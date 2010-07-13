/*
 * BranchSegment.h
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */

#ifndef BRANCHSEGMENT_COPPER_H_
#define BRANCHSEGMENT_COPPER_H_

#include <vector>
using namespace std;

#include "RateModel.h"

#include "vector_node_object.h"

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

class BranchSegment{
	private:
		double duration;
		int period;
		RateModel * model;
		vector<int> fossilareaindices;
		int startdistint;

	public:
		BranchSegment(double dur,int per);
		void setModel(RateModel * mod);
		//void setStartDist(vector<int> sd);
		void clearStartDist();
		double getDuration();
		int getPeriod();
		//vector<int> getStartDist();
		void set_start_dist_int(int d);
		int get_start_dist_int();
		RateModel * getModel();
		vector<int> getFossilAreas();
		void setFossilArea(int area);
#ifdef BIGTREE
		VectorNodeObject<mpfr_class> * distconds;
		VectorNodeObject<mpfr_class> alphas;
		VectorNodeObject<mpfr_class> * ancdistconds;//for ancestral state reconstructions
#else
		VectorNodeObject<double> * distconds;
		VectorNodeObject<double> alphas; // alpha for the entire branch -- stored in the 0th segment for anc calc
		VectorNodeObject<double> seg_sp_alphas; // alpha for this specific segment, stored for the stoch map
		VectorNodeObject<double> seg_sp_stoch_map_revB_time; //segment specific rev B, combining the tempA and the ENLT
		VectorNodeObject<double> seg_sp_stoch_map_revB_number; //segment specific rev B, combining the tempA and the ENLT
		VectorNodeObject<double> * ancdistconds;//for ancestral state reconstructions
#endif
};

#endif /* BRANCHSEGMENT_H_ */
