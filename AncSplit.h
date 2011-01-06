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

#include "superdouble.h"
#include "RateModel.h"

class AncSplit{
	private:
		RateModel * model;
		double weight;
		Superdouble likelihood;

	public:
		AncSplit(RateModel * mod,int,int,int,Superdouble);
		RateModel * getModel();
		double getWeight();
		Superdouble getLikelihood();
		void setLikelihood(Superdouble li);
		int ancdistint;
		int ldescdistint;
		int rdescdistint;
};


#endif /* ANCSPLIT_H_ */
