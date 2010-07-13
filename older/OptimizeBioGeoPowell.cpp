/*
 * OptimizeBioGeoPowell.cpp
 *
 *  Created on: Aug 21, 2009
 *      Author: smitty
 */

#include "OptimizeBioGeoPowell.h"
#include "BioGeoTree.h"
#include "RateModel.h"

#include <math.h>
#include <stdio.h>
#include <float.h>
#include <limits>
using namespace std;

namespace {
	inline void shft3(double &a, double &b, double &c, const double d)
	{
		a=b;
		b=c;
		c=d;
	}
	inline const double SQR(const double a) {return a*a;}
	inline double MAX(const double &a, const double &b)
	        {return b > a ? (b) : double(a);}
	inline void SWAP(double &a, double &b)
		{double dum=a; a=b; b=dum;}
}

OptimizeBioGeoPowell::OptimizeBioGeoPowell(BioGeoTree * intree, RateModel * inrm,bool marg){
	NDIM=2;
	FTOL=0.001;//stopping criterion fraction
	ITMAX = 50;
	CGOLD = 0.3819660;
	GOLD = 1.618034;
	ZEPS = 0.000000001;
	TINY = 1.0E-20;
	GLIMIT = 100.0;
	rm = inrm;
	tree = intree;
}

vector<double> OptimizeBioGeoPowell::optimize_global_dispersal_extinction(){
	vector<double> p_d(2,0.01);
	vector<double> p(NDIM,0.01);
	vector<vector<double> > xi = vector<vector<double> > (NDIM,vector<double>(NDIM));
	for (int i=0;i<NDIM;i++)
		for (int j=0;j<NDIM;j++)
			xi[i][j]=(i == j ? 1.0 : -1.0);
	powell(p,xi,FTOL,iter,fret);
/*	cout << "Iterations: " << iter << endl << endl;;
	cout << "Minimum found at: " << endl;
	cout << fixed << setprecision(6);
	for (int i=0;i<NDIM;i++) cout << setw(12) << p[i];
	cout << endl << endl << "Minimum function value = ";
	cout << setw(12) << fret << endl << endl;
	cout << "True minimum of function is at:" << endl;
	cout << setw(12) << 1.0 << setw(12) << 2.0 << setw(12) << 3.0 << endl;*/
	return p;
}

void OptimizeBioGeoPowell::powell(vector<double> &p, vector<vector<double> > &xi, const double ftol, int &iter,double &fret){
	ITMAX=200;
	TINY=1.0e-25;
	int i,j,ibig;
	double del,fp,fptt,t;

	int n=p.size();
	vector<double> pt(n),ptt(n),xit(n);
	fret=func(p);
	for (j=0;j<n;j++) pt[j]=p[j];
	for (iter=0;;++iter) {
		fp=fret;
		ibig=0;
		del=0.0;
		for (i=0;i<n;i++) {
			for (j=0;j<n;j++) xit[j]=xi[j][i];
			fptt=fret;
			linmin(p,xit,fret);
			if (fptt-fret > del) {
				del=fptt-fret;
				ibig=i+1;
			}
		}
		if (2.0*(fp-fret) <= ftol*(fabs(fp)+fabs(fret))+TINY) {
			return;
		}
		if (iter == ITMAX) cout << "powell exceeding maximum iterations." << endl;
		for (j=0;j<n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=func(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin(p,xit,fret);
				for (j=0;j<n;j++) {
					xi[j][ibig-1]=xi[j][n-1];
					xi[j][n-1]=xit[j];
				}
			}
		}
	}
}

void OptimizeBioGeoPowell::linmin(vector<double> &p, vector<double> &xi, double &fret){
	int j;
	const double TOL=1.0e-8;
	double xx,xmin,fx,fb,fa,bx,ax;

	int n=p.size();
	ncom=n;
	pcom_p=new vector<double>(n);
	xicom_p=new vector<double>(n);
	//nrfunc=func;
	vector<double> &pcom=*pcom_p,&xicom=*xicom_p;
	for (j=0;j<n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(ax,xx,bx,fa,fx,fb);
	fret=brent(ax,xx,bx,TOL,xmin);
	for (j=0;j<n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	delete xicom_p;
	delete pcom_p;
}

void OptimizeBioGeoPowell::mnbrak(double &ax, double &bx, double &cx, double &fa, double &fb, double &fc){
	const double GOLD=1.618034,GLIMIT=100.0,TINY=1.0e-20;
	double ulim,u,r,q,fu;

	fa=f1dim(ax);
	fb=f1dim(bx);
	if (fb > fa) {
		SWAP(ax,bx);
		SWAP(fb,fa);
	}
	cx=bx+GOLD*(bx-ax);
	fc=f1dim(cx);
	while (fb > fc) {
		r=(bx-ax)*(fb-fc);
		q=(bx-cx)*(fb-fa);
		u=bx-((bx-cx)*q-(bx-ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=bx+GLIMIT*(cx-bx);
		if ((bx-u)*(u-cx) > 0.0) {
			fu=f1dim(u);
			if (fu < fc) {
				ax=bx;
				bx=u;
				fa=fb;
				fb=fu;
				return;
			} else if (fu > fb) {
				cx=u;
				fc=fu;
				return;
			}
			u=cx+GOLD*(cx-bx);
			fu=f1dim(u);
		} else if ((cx-u)*(u-ulim) > 0.0) {
			fu=f1dim(u);
			if (fu < fc) {
				shft3(bx,cx,u,cx+GOLD*(cx-bx));
				shft3(fb,fc,fu,f1dim(u));
			}
		} else if ((u-ulim)*(ulim-cx) >= 0.0) {
			u=ulim;
			fu=f1dim(u);
		} else {
			u=cx+GOLD*(cx-bx);
			fu=f1dim(u);
		}
		shft3(ax,bx,cx,u);
		shft3(fa,fb,fc,fu);
	}
}

double OptimizeBioGeoPowell::brent(const double ax, const double bx, const double cx,
	const double tol, double &xmin)
{
	const int ITMAX=100;
	const double CGOLD=0.3819660;
	const double ZEPS=numeric_limits<double>::epsilon()*1.0e-3;
	int iter;
	double a,b,d=0.0,etemp,fu,fv,fw,fx;
	double p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;
	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=f1dim(x);
	for (iter=0;iter<ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=f1dim(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			shft3(v,w,x,u);
			shft3(fv,fw,fx,fu);
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	cout << "Too many iterations in brent" <<endl;;
	xmin=x;
	return fx;
}

double OptimizeBioGeoPowell::f1dim(const double x)
{
	int j;
	vector<double> xt(ncom);
	vector<double> &pcom=*pcom_p,&xicom=*xicom_p;
	for (j=0;j<ncom;j++)
		xt[j]=pcom[j]+x*xicom[j];
	return func(xt);
}

double OptimizeBioGeoPowell::func(vector<double> &p){
	double like;
	if (p[0] <= 0 || p[1] <= 0){
		like = 10000000;
	}else{
		rm->setup_D(p[0]);
		rm->setup_E(p[1]);
		rm->setup_Q();
		tree->update_default_model(rm);
		like = -log(tree->eval_likelihood(true));
		if(like == std::numeric_limits<double>::infinity())
			like = 10000000;
	}
	//cout <<p[0] << " " << p[1] << " " <<  like << endl;
	return like;
}

double OptimizeBioGeoPowell::SIGN(double a,double b){
	double tempS;
	if(b>=0.0)
		tempS=fabs(a);
	else
		tempS=-fabs(a);
	return tempS;
}



