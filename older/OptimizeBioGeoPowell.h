/*
 * OptimizeBioGeoPowell.h
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#ifndef OPTIMIZEBIOGEOPOWELL_H_
#define OPTIMIZEBIOGEOPOWELL_H_

#include <vector>
using namespace std;

#include "BioGeoTree.h"
#include "RateModel.h"

class OptimizeBioGeoPowell{
	private:
		BioGeoTree * tree;
		RateModel * rm;
		int ITMAX;
		double CGOLD;
		double GOLD;
		double ZEPS;
		double TINY;
		double GLIMIT;
		double fret;
		int NDIM;
		double FTOL;
		int iter;
		int ncom;
		vector<double> * pcom_p;
		vector<double> * xicom_p;
		double xmin;
		double * mainP;
		double mainFret;
		double ax;
		double bx;
		double cx;
		//from mnbrak
		double mnfa;
		double mnfb;
		double mnfc;
		void powell(vector<double> &p, vector<vector<double> > &xi, const double ftol, int &iter,double &fret);
		void linmin(vector<double> &p, vector<double> &xi, double &fret);
		double f1dim(const double x);
		double brent(const double ax, const double bx, const double cx,
				const double tol, double &xmin);
		void mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc);
		double func(vector<double> &p);
		double SIGN(double a,double b);

	public:
		OptimizeBioGeoPowell(BioGeoTree * intree,RateModel * inrm, bool marg);
		vector<double> optimize_global_dispersal_extinction();


};
#endif /* OPTIMIZEBIOGEOPOWELL_H_ */
