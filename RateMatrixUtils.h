/*
 * RateMatrixUtils.h
 *
 *  Created on: Aug 13, 2009
 *      Author: Stephen A. Smith
 */

/*
  This is essentially a utility file that stores matrix utility functions dealing
  with the rate matrix.
 */

#ifndef RATEMATRIXUTILS_H_
#define RATEMATRIXUTILS_H_
#include "RateModel.h"
#include "AncSplit.h"
#include <vector>
using namespace std;

/*
  most of these utilities calculate simple math on matrices and vectors
  and they should only be used if there is the chance that there will 
  be a null vector or matrix. otherwise, c++ numerics library should
  be used for speed.
 */
double calculate_vector_double_sum(vector<double> & in);
int calculate_vector_int_sum(vector<int> * in);
int get_vector_int_index_from_multi_vector_int(vector<int> * in, vector<vector<int> > * in2);

/*
  used for creating and enumerate distributions
 */
int calculate_vector_int_sum_xor(vector<int> & in,vector<int> & in2);
int locate_vector_int_single_xor(vector<int> & in, vector<int> & in2);

/*
  given the rate model and the distribution, this should return all the 
  ancestral splits and should ONLY be used for calculating ancestral 
  state values. otherwise the one returning only the ints should be returned
 */
vector<AncSplit> iter_ancsplits(RateModel *rm, vector<int> & dist);

/*
  like the above function but without using AncSplits object, and should be
  used for everything but ancestral state reconstruction.
  output is the leftdists and rightdists which are the index of the distribution
  in the ratemodel->getdists
 */
void iter_ancsplits_just_int(RateModel *rm, vector<int> & dist,vector<int> & leftdists, vector<int> & rightdists, double & weight);

/*
  simple printing functions
 */
void print_vector_int(vector<int> & in);
void print_vector_double(vector<double> & in);

/*
  used for generating all the distributions with maximium number of areas involved
  would be designated in the config file
 */
vector<vector<int> > generate_dists_from_num_max_areas(int totalnum,int numareas);

/*
  used for processing custom rate matrix config files designated in the main config file
 */
vector<vector<vector<double> > > processRateMatrixConfigFile(string filename, int numareas, int nperiods);

/*
  all of these are for sparse matrix calculation and are used for the fortran methods

  WILL PROBABLY CHANGE WHEN MOVED TO C++ CLASSES FOR MATEXP
 */
int get_size_for_coo(vector<vector<double> > & inm, double t);
void convert_matrix_to_coo_for_fortran(vector<vector<double> > & inmatrix, double t, int * ia, int * ja, double * a);
void convert_matrix_to_coo_for_fortran_vector(vector<vector<double> > & inmatrix, vector<int> & ia, vector<int> & ja, vector<double> & a);
void convert_matrix_to_single_row_for_fortran(vector<vector<double> > & inmatrix, double t, double * H);
vector<int> get_columns_for_sparse(vector<double> &,RateModel *);

/*
	this is for pthread sparse columns

	PRETTY MUCH USELESS UNTIL IMPLEMENTED MATEXP STUFF IN C++
 */
struct sparse_thread_data{
	int thread_id;
	vector<int> columns;
	vector<vector<double> > presults;
	RateModel * rm;
	double t;
	int period;
};

void * sparse_column_pmatrix_pthread_go(void *threadarg);

/*
 * REQUIRES BOOST -- UNCOMMENT TO REACTIVATE
 */
/*
 * this constructs a P matrix using a Q matrix
 */

//vector< vector<double> > QMatrixToPmatrix(vector< vector<double> > & Q, double t);
//void calcMatExp(int * ia,int * ja, double * a, int n);


#endif /* RATEMATRIXUTILS_H_ */
