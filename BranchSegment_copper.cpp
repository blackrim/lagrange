/*
 * BranchSegment.cpp
 *
 *  Created on: Aug 16, 2009
 *      Author: smitty
 */


#include "BranchSegment_copper.h"
#include "RateModel.h"

#include <vector>
using namespace std;

#include "vector_node_object.h"

BranchSegment::BranchSegment(double dur,int per){
	duration = dur;
	period = per;
	startdistint = -666;
}

void BranchSegment::setModel(RateModel * mod){
	model = mod;
}

/*void BranchSegment::setStartDist(vector<int> sd){
	startdist = sd;
}*/

void BranchSegment::clearStartDist(){
	//startdist.clear();
	startdistint = -666; //null is -666
}

double BranchSegment::getDuration(){
	return duration;
}

int BranchSegment::getPeriod(){
	return period;
}
/*
vector<int> BranchSegment::getStartDist(){
	return startdist;
}*/

void BranchSegment::set_start_dist_int(int d){
	startdistint = d;
}

int BranchSegment::get_start_dist_int(){
	return startdistint;
}

RateModel * BranchSegment::getModel(){
	return model;
}

vector<int> BranchSegment::getFossilAreas(){
	return fossilareaindices;
}

void BranchSegment::setFossilArea(int area){
	fossilareaindices.push_back(area);
}
