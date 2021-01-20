/*
 * RateMatrix.cpp
 *
 *  Created on: Aug 14, 2009
 *      Author: smitty
 */

#define VERBOSE false

#include "RateModel.h"
#include "RateMatrixUtils.h"
#include "Utils.h"
//#include "AncSplit.h"

#include <pthread.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <map>
#include <math.h>
using namespace std;

#include <armadillo>
using namespace arma;

RateModel::RateModel(int na, bool ge, vector<double> pers, bool sp):
	globalext(ge),nareas(na),numthreads(0),periods(pers),sparse(sp){}

void RateModel::set_nthreads(int nthreads){
	numthreads = nthreads;
}

int RateModel::get_nthreads(){
	return numthreads;
}

void RateModel::setup_dists(){
	map< int, vector<int> > a = iterate_all_bv(nareas);
	if (globalext){
		vector<int> empt;
		for (unsigned int i=0;i<a[0].size();i++){
			empt.push_back(0);
		}
		dists.push_back(empt);
	}
	map<int, vector<int> >::iterator pos;
	for (pos = a.begin(); pos != a.end(); ++pos){
		int f = pos->first;
		dists.push_back(a[f]);
	}
	/*
	 calculate the distribution map
	 */
	for(unsigned int i=0;i<dists.size();i++){
		distsintmap[dists[i]] = i;
	}
	for(unsigned int i=0;i<dists.size();i++){
		intdistsmap[i] = dists[i];
	}
	/*
	 * precalculate the iterdists
	 */
	iter_all_dist_splits();

	/*
	 print out a visual representation of the matrix
	 */
	if (VERBOSE){
		cout << "dists" <<endl;
		for (unsigned int j=0; j< dists.size(); j++){
			cout << j << " ";
			for (unsigned int i=0;i<dists[j].size();i++){
				cout << dists[j][i];
			}
			cout << endl;
		}
	}
}

/*
 * need to make a generator function for setting distributions
 */
void RateModel::setup_dists(vector<vector<int> > indists, bool include){
	if(include == true){
		dists = indists;
//		if(calculate_vector_int_sum(&dists[0]) > 0){
		if(accumulate(dists[0].begin(),dists[0].end(),0) > 0){
			vector<int> empt;
			for (unsigned int i=0;i<dists[0].size();i++){
				empt.push_back(0);
			}
			dists.push_back(empt);
		}
	}else{//exclude is sent
		vector<int> empt;
		for (int i=0;i<nareas;i++){
			empt.push_back(0);
		}
		dists.push_back(empt);

		map< int, vector<int> > a = iterate_all_bv(nareas);
		map<int, vector<int> >::iterator pos;
		for (pos = a.begin(); pos != a.end(); ++pos){
			int f = pos->first;
			bool inh = false;
			for(unsigned int j=0;j<indists.size();j++){
				//if(is_vector_int_equal_to_vector_int(indists[j],a[f])){
				if(indists[j]==a[f]){
					inh = true;
				}
			}
			if(inh == false)
				dists.push_back(a[f]);
		}
	}
	/*
	 calculate the distribution map
	 */
	for(unsigned int i=0;i<dists.size();i++){
		distsintmap[dists[i]] = i;
	}
	for(unsigned int i=0;i<dists.size();i++){
		intdistsmap[i] = dists[i];
	}
	/*
	precalculate the iterdists
	 */
	iter_all_dist_splits();

	/*
	 print out a visual representation of the matrix
	 */
	if (VERBOSE){
		cout << "dists" <<endl;
		for (unsigned int j=0; j< dists.size(); j++){
			cout << j << " ";
			for (unsigned int i=0;i<dists[j].size();i++){
				cout << dists[j][i];
			}
			cout << endl;
		}
	}
}


/*
void RateModel::remove_dist(vector<int> dist);*/

/*
 * just give Dmask a bunch of ones
 * specify particular ones in the Dmask_cell
 */
void RateModel::setup_Dmask(){
	vector<double> cols(nareas, 1);
	vector< vector<double> > rows(nareas, cols);
	Dmask = vector< vector< vector<double> > > (periods.size(), rows);
}

void RateModel::set_Dmask_cell(int period, int area, int area2, double prob, bool sym){
	Dmask[period][area][area2] = prob;
	if (sym)
		Dmask[period][area2][area] = prob;
}

void RateModel::setup_D(double d){
	vector<double> cols(nareas, 1*d);
	vector< vector<double> > rows(nareas, cols);
	D = vector< vector< vector<double> > > (periods.size(), rows);
	for (unsigned int i=0;i<D.size();i++){
		for (unsigned int j=0;j<D[i].size();j++){
			D[i][j][j] = 0.0;
			for (unsigned int k=0;k<D[i][j].size();k++){
				D[i][j][k] = D[i][j][k] * Dmask[i][j][k];
			}
		}
	}
	if (VERBOSE){
		cout << "D" <<endl;
		for (unsigned int i=0;i<D.size();i++){
			for (unsigned int j=0;j<D[i].size();j++){
				for (unsigned int k=0;k<D[i][j].size();k++){
					cout << D[i][j][k] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
	}
}

/*
 * this is for estimating the D matrix
 * in this case, a setup dmatrix is being sent and
 * the dmask is then applied to it
 */

void RateModel::setup_D_provided(double d, vector< vector< vector<double> > > & D_mask_in){
	vector<double> cols(nareas, 1*d);
	vector< vector<double> > rows(nareas, cols);
	D = vector< vector< vector<double> > > (periods.size(), rows);
	for (unsigned int i=0;i<D.size();i++){
		for (unsigned int j=0;j<D[i].size();j++){
			D[i][j][j] = 0.0;
			for (unsigned int k=0;k<D[i][j].size();k++){
				D[i][j][k] = D[i][j][k] * Dmask[i][j][k]*D_mask_in[i][j][k];
			}
		}
	}
	if (VERBOSE){
		cout << "D" <<endl;
		for (unsigned int i=0;i<D.size();i++){
			for (unsigned int j=0;j<D[i].size();j++){
				for (unsigned int k=0;k<D[i][j].size();k++){
					cout << D[i][j][k] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
	}
}

void RateModel::setup_E(double e){
	vector<double> cols(nareas, 1*e);
	E = vector<vector<double> > (periods.size(), cols);
}

void RateModel::set_Qdiag(int period){
	for (unsigned int i=0;i<dists.size();i++){
		double sum =(calculate_vector_double_sum(Q[period][i]) - Q[period][i][i]) * -1.0;
		Q[period][i][i] = sum;
	}
}

void RateModel::setup_Q(){
	vector<double> cols(dists.size(), 0);
	vector< vector<double> > rows(dists.size(), cols);
	Q = vector< vector< vector<double> > > (periods.size(), rows);
	for(unsigned int p=0; p < Q.size(); p++){//periods
		for(unsigned int i=0;i<dists.size();i++){//dists
			//int s1 = calculate_vector_int_sum(&dists[i]);
			int s1 = accumulate(dists[i].begin(),dists[i].end(),0);
			if(s1 > 0){
				for(unsigned int j=0;j<dists.size();j++){//dists
					int sxor = calculate_vector_int_sum_xor(dists[i], dists[j]);
					if (sxor == 1){
						//int s2 = calculate_vector_int_sum(&dists[j]);
						int s2 = accumulate(dists[j].begin(),dists[j].end(),0);
						int dest = locate_vector_int_single_xor(dists[i],dists[j]);
						double rate = 0.0;
						if (s1 < s2){
							for (unsigned int src=0;src<dists[i].size();src++){
								if(dists[i][src] != 0){
									rate += D[p][src][dest] ;//* Dmask[p][src][dest];
								}
							}
						}else{
							rate = E[p][dest];
						}
						Q[p][i][j] = rate;
					}
				}
			}
		}
		set_Qdiag(p);
	}
	/*
	 * sparse needs to be transposed for matrix exponential calculation
	 */
	if(sparse == true){
		vector<double> cols(dists.size(), 0);
		vector< vector<double> > rows(dists.size(), cols);
		QT = vector< vector< vector<double> > > (periods.size(), rows);
		for(unsigned int p=0; p < QT.size(); p++){//periods
			for(unsigned int i=0;i<dists.size();i++){//dists
				for(unsigned int j=0;j<dists.size();j++){//dists
					QT[p][j][i] = Q[p][i][j];
				}
			}
		}
		//setting up the coo numbs
		nzs = vector<int>(Q.size(),0);
		for(unsigned int p=0; p < Q.size(); p++){//periods
			nzs[p] = get_size_for_coo(Q[p],1);
		}
		//setup matrix
		ia_s.clear();
		ja_s.clear();
		a_s.clear();
		for(unsigned int p=0; p < Q.size(); p++){//periods
			vector<int> ia = vector<int>(nzs[p]);
			vector<int> ja = vector<int>(nzs[p]);
			vector<double> a = vector<double>(nzs[p]);
			convert_matrix_to_coo_for_fortran_vector(QT[p],ia,ja,a);//need to multiply these all these by t
			ia_s.push_back(ia);
			ja_s.push_back(ja);
			a_s.push_back(a);
		}
	}
	if(VERBOSE){
		cout << "Q" <<endl;
		for (unsigned int i=0;i<Q.size();i++){
			for (unsigned int j=0;j<Q[i].size();j++){
				for (unsigned int k=0;k<Q[i][j].size();k++){
					cout << Q[i][j][k] << " ";
				}
				cout << endl;
			}
			cout << endl;
		}
	}
}

vector<vector<double > > RateModel::setup_arma_P(int period, double t, bool store_p_matrices){
	int m = Q[period].size();
	mat A(m, m);
	for(int i=0;i<m;i++){
		for(int j=0;j<m;j++){
			A.at(i,j) = Q[period][i][j];
		}
	}
	mat B = expmat(A);
	vector<vector<double> > p (Q[period].size(), vector<double>(Q[period].size()));
	for(int i=0;i<m;i++){
		for(int j=0;j<m;j++){
			p[i][j] = B.at(i,j);
		}
	}
	if(store_p_matrices == true){
		stored_p_matrices[period][t] = p;
	}
	return p;
}

vector<vector<vector<int> > > RateModel::iter_dist_splits(vector<int> & dist){
	vector< vector <vector<int> > > ret;
	vector< vector<int> > left;
	vector< vector<int> > right;
	if(accumulate(dist.begin(),dist.end(),0) == 1){
		left.push_back(dist);
		right.push_back(dist);
	}
	else{
		for(unsigned int i=0;i<dist.size();i++){
			if (dist[i]==1){
				vector<int> x(dist.size(),0);
				x[i] = 1;
				int cou = count(dists.begin(),dists.end(),x);
				if(cou > 0){
					left.push_back(x);right.push_back(dist);
					left.push_back(dist);right.push_back(x);
					vector<int> y;
					for(unsigned int j=0;j<dist.size();j++){
						if(dist[j]==x[j]){
							y.push_back(0);
						}else{
							y.push_back(1);
						}
					}
					int cou2 = count(dists.begin(),dists.end(),y);
					if(cou2 > 0){
						left.push_back(x);right.push_back(y);
						if(accumulate(y.begin(),y.end(),0) > 1){
							left.push_back(y);right.push_back(x);
						}
					}
				}
			}
		}
	}
	if(VERBOSE){
		cout << "LEFT" << endl;
		for(unsigned int i = 0; i< left.size() ; i++ ){
			print_vector_int(left[i]);
		}
		cout << "RIGHT" << endl;
		for(unsigned int i = 0; i< right.size() ; i++ ){
			print_vector_int(right[i]);
		}
	}
	ret.push_back(left);
	ret.push_back(right);
	return ret;
}

void RateModel::iter_all_dist_splits(){
	for(unsigned int i=0;i<dists.size();i++){
		iter_dists[dists[i]] = iter_dist_splits(dists[i]);
	}
}


vector<vector<int> > * RateModel::getDists(){
	return &dists;
}

map<vector<int>,int> * RateModel::get_dists_int_map(){
	return &distsintmap;
}

map<int,vector<int> > * RateModel::get_int_dists_map(){
	return &intdistsmap;
}

vector<vector<vector<int> > > * RateModel::get_iter_dist_splits(vector<int> & dist){
	return &iter_dists[dist];
}

int RateModel::get_num_areas(){return nareas;}

int RateModel::get_num_periods(){return periods.size();}

vector< vector< vector<double> > > & RateModel::get_Q(){
	return Q;
}

inline int signof(double d)
{
   return d >= 0 ? 1 : -1;
}

inline double roundto(double in){
	return floor(in*(1000)+0.5)/(1000);
}

/*
 * this should be used to caluculate the eigenvalues and eigenvectors
 * as U * Q * U-1 -- eigen decomposition
 *
 * this should use the armadillo library
 */
bool RateModel::get_eigenvec_eigenval_from_Q(cx_mat * eigval, cx_mat * eigvec, int period){
	mat tQ(int(Q[period].size()),int(Q[period].size())); tQ.fill(0);
	for(unsigned int i=0;i<Q[period].size();i++){
		for(unsigned int j=0;j<Q[period].size();j++){
			tQ(i,j) = Q[period][i][j];
			//cout << Q[0][i][j] << " ";
		}
		//cout << endl;
	}
	//cout << endl;
	cx_colvec eigva;
	cx_mat eigve;
	eig_gen(eigva,eigve,tQ);
	bool isImag = false;
	for(unsigned int i=0;i<Q[period].size();i++){
		for(unsigned int j=0;j<Q[period].size();j++){
			if(i==j)
				(*eigval)(i,j) = eigva(i);
			else
				(*eigval)(i,j) = 0;
			(*eigvec)(i,j) = eigve(i,j);
			if(imag((*eigvec)(i,j)) > 0 || imag((*eigval)(i,j)))
				isImag = true;
		}
	}
	if (VERBOSE) {
        cout << eigva << endl;
        cout << tQ - ((*eigvec) * (*eigval) * inv(*eigvec)) << endl;
    }
	return isImag;
}

//trying not to use octave at the moment
/*
 * 
 * bool RateModel::get_eigenvec_eigenval_from_Q_octave(ComplexMatrix * eigval, ComplexMatrix * eigvec, int period){
	ComplexMatrix tQ = ComplexMatrix(int(Q[period].size()),int(Q[period].size()));
	for(unsigned int i=0;i<Q[period].size();i++){
		for(unsigned int j=0;j<Q[period].size();j++){
			tQ(i,j) = Q[period][i][j];
	//		cout << Q[0][i][j] << " ";
		}
	//	cout << endl;
	}
	//cout << endl;
	EIG eig = EIG(tQ);
	bool isImag = false;
	for(unsigned int i=0;i<Q[period].size();i++){
		for(unsigned int j=0;j<Q[period].size();j++){
			if(i==j){
				(*eigval)(i,j) = eig.eigenvalues()(i);
			}else{
				(*eigval)(i,j) = 0;
			}
			(*eigvec)(i,j) = eig.eigenvectors()(i,j);
			if(imag((*eigvec)(i,j)) > 0 || imag((*eigval)(i,j)))
				isImag = true;
		}
	}
	return isImag;
	//cout <<(eig.eigenvalues() * eig.eigenvectors()) << endl;
}*/

/**/


//REQUIRES BOOST AND IS SLOWER BUT TO ACTIVATE UNCOMMENT
/*vector<vector<double > > RateModel::setup_P(int period, double t){
	//
	//return P, the matrix of dist-to-dist transition probabilities,
	//from the model's rate matrix (Q) over a time duration (t)
	//
	vector<vector<double> > p = QMatrixToPmatrix(Q[period], t);

	//filter out impossible dists
	//vector<vector<int> > dis = enumerate_dists();
	for (unsigned int i=0;i<dists.size();i++){
		//if (calculate_vector_int_sum(&dists[i]) > 0){
		if(accumulate(dists[i].begin(),dists[i].end(),0) > 0){
			for(unsigned int j=0;j<dists[i].size();j++){
				if(dists[i][j]==1){//present
					double sum1 =calculate_vector_double_sum(Dmask[period][j]);
					double sum2 = 0.0;
					for(unsigned int k=0;k<Dmask[period].size();k++){
						sum2 += Dmask[period][k][j];
					}
					if(sum1+sum2 == 0){
						for(unsigned int k=0;k<p[period].size();k++){
							p[period][k] = p[period][k]*0.0;
						}
						break;
					}
				}
			}
		}
	}
	if(VERBOSE){
		cout << "p " << period << " "<< t << endl;
		for (unsigned int i=0;i<p.size();i++){
			for (unsigned int j=0;j<p[i].size();j++){
				cout << p[i][j] << " ";
			}
			cout << endl;
		}
	}
	return p;
}*/



