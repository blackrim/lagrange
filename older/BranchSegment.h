/*
 * BranchSegment.h
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */

#ifndef BRANCHSEGMENT_H_
#define BRANCHSEGMENT_H_

#include <vector>
using namespace std;

#include "RateModel.h"

#include <Utils/BppVector.h>
using namespace bpp;

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
		Vector<double> * distconds;
		Vector<double> alphas;
		Vector<double> * ancdistconds;//for ancestral state reconstructions
};

#endif /* BRANCHSEGMENT_H_ */
