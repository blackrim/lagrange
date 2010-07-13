/*
 *  BayesianBioGeo.cpp
 *  lagrange_cpp
 *
 *  Created by Stephen Smith on 1/13/10.
 *  Copyright 2010 Yale University. All rights reserved.
 *
 */

#include "BayesianBioGeo.h"

#include <fstream>
#include <map>
#include <vector>
#include <math.h>
#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
using namespace std;

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "RateModel.h"
#include "BioGeoTree_copper.h"

namespace {
	inline double MIN(const double &a, const double &b)
	{return b < a ? (b) : double(a);}
}

BayesianBioGeo::BayesianBioGeo(BioGeoTree_copper * intree,RateModel * inrm, bool marg, int gen):
	tree(intree),rm(inrm), gens(gen), marginal(marg){
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
}

void BayesianBioGeo::run_global_dispersal_extinction(){
	double prevlike = 0;
	double prevprior = 0;
	double prevpost = 0;
	double curlike = 0;
	double curpost = 0;
	double curprior = 0;

	vector<double> sliding(2);
	vector<double> trials(2);
	vector<double> success(2);
	sliding[0] = 0.5;sliding[1] = 0.5;
	for(unsigned int i=0;i<trials.size();i++){trials[i] = 0;}
	for(unsigned int i=0;i<success.size();i++){success[i] = 0;}
	
	int rot = 0;
	double hastings = 1;
	
	ofstream outfile ("test.txt");

	params = vector<double>(2);
	prevparams = vector<double>(2);
	for(unsigned int i=0;i<params.size();i++){params[i] = 0.1;}
	rm->setup_D(0.1);
	rm->setup_E(0.1);
	rm->setup_Q();
	tree->update_default_model(rm);
	prevlike = -tree->eval_likelihood(marginal);
	prevprior = 1;
	for(unsigned int i=0;i<params.size();i++){prevprior *= calculate_pdf(params[i]);}
	cout << "intial like: " << prevlike << endl;
	int iter = 0;
	while(iter < gens){
		double dispersal=params[0];
		double extinction=params[1];
		rm->setup_D(dispersal);
		rm->setup_E(extinction);
		rm->setup_Q();
		tree->update_default_model(rm);
		curlike = -tree->eval_likelihood(marginal);
		/*
		 * calcprior
		 */
		curprior = 1;
		for(unsigned int i=0;i<params.size();i++){curprior *= calculate_pdf(params[i]);}
		
		/*
		 * check to keep
		 */
		double testr = gsl_ran_flat(r,0,1);
		double test = MIN(1,hastings * (curprior/prevprior) * (exp(curlike-prevlike)));
		
		if (iter > 100)
			trials[rot] += 1;
		//cout << "- "<< prevlike << " " <<  curlike  << " " << hastings << " " << test << " " << testr;
		//cout << " " << rot << " " << dispersal << " " << extinction << endl;
		if (testr < test){
			prevprior = curprior;
			prevlike = curlike;
			prevpost = curpost;
			prevparams[rot] = params[rot];
			if (iter > 100)
				success[rot] += 1;
		}
		/*
		 * pick next params
		 */
		//params[rot] = calculate_sliding(prevparams[rot],sliding[rot]);
		params[rot] = calculate_sliding_log(prevparams[rot],sliding[rot], &hastings);
		if(iter%5 == 0){
			rot += 1;
		}
		if (rot == 2){
			rot = 0;
		}
		if(iter%100 == 0 && iter > 100){
			cout << iter << " " << prevlike;
			for(unsigned int i=0;i<params.size();i++){cout << " " << prevparams[i];}
			for(unsigned int i=0;i<params.size();i++){cout << " " << success[i]/trials[i];}
			cout << endl;
			outfile << iter << "\t" << prevlike;
			for(unsigned int i=0;i<params.size();i++){outfile << "\t" << prevparams[i];}
			outfile << endl;
		}
		iter++;
	}
	outfile.close();
}

double BayesianBioGeo::calculate_pdf(double value){
	return gsl_ran_flat_pdf (value, 0,100.);
}

double BayesianBioGeo::calculate_sliding(double value, double sliding){
	//return abs(gsl_ran_flat(r,value - (sliding/2),value + (sliding/2)));
	return abs(value + gsl_ran_gaussian(r,sliding));
	//cauchy
}

double BayesianBioGeo::calculate_sliding_log(double value, double sliding, double * hastings){
	double newv = log(value) + gsl_ran_gaussian(r,sliding);
	newv = abs(exp(newv));
	*hastings = 1 * newv/value;
	return newv;
}
