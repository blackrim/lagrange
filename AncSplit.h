/*
 * AncSplit.h
 *
 *  Created on: Aug 15, 2009
 *      Author: smitty
 */

/*
  The AncSplit (ancestral split) class is used to store the likelihoods and 
  distributions associated with an ancstral split. 

  This should only be used for ancestral state calculation as there is no
  need to store the likelihoods for each state (other than distconds) when 
  calculating the likeklihood.  
 */

#ifndef ANCSPLIT_H_
#define ANCSPLIT_H_

#include <vector>
//using namespace std;

#include "RateModel.h"

#ifdef BIGTREE
#include "gmpfrxx/gmpfrxx.h"
#endif

class AncSplit{
	private:
		RateModel * model;
		double weight;
#ifdef BIGTREE
		mpfr_class likelihood;
#else
		double likelihood;
#endif

	public:
		AncSplit(RateModel * mod,int,int,int,double);
		RateModel * getModel();
		double getWeight();
#ifdef BIGTREE
		mpfr_class getLikelihood();
		void setLikelihood(mpfr_class li);
#else
		double getLikelihood();
		void setLikelihood(double li);
#endif
		int ancdistint;
		int ldescdistint;
		int rdescdistint;
};


#endif /* ANCSPLIT_H_ */
