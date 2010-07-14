/*
 * AncSplit.cpp
 *
 *  Created on: Aug 15, 2009
 *      Author: Stephen A. Smith
 */

/*
  The AncSplit (ancestral split) class is used to store the likelihoods and 
  distributions associated with an ancstral split. 

  This should only be used for ancestral state calculation as there is no
  need to store the likelihoods for each state (other than distconds) when 
  calculating the likeklihood.  
 */

#include "AncSplit.h"
#include "RateModel.h"

#include <vector>
using namespace std;

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif


AncSplit::AncSplit(RateModel * mod,int dist,int ldesc,int rdesc,double we):model(mod),weight(we),
		likelihood(0),ancdistint(dist),ldescdistint(ldesc),rdescdistint(rdesc){}

RateModel * AncSplit::getModel(){
	return model;
}

double AncSplit::getWeight(){
	return weight;
}


#ifdef BIGTREE
mpfr_class AncSplit::getLikelihood(){
	return likelihood;
}

void AncSplit::setLikelihood(mpfr_class li){
	likelihood = li;
}
#else
double AncSplit::getLikelihood(){
	return likelihood;
}

void AncSplit::setLikelihood(double li){
	likelihood = li;
}
#endif
