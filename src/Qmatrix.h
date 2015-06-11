/*
 * Qmatrix.h
 *
 *  Created on: Nov 8, 2011
 *      Author: Tuomas Rajala <tuomas.rajala@jyu.fi>
 */

#ifndef QMATRIX_H_
#define QMATRIX_H_
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <vector>

class Qmatrix {
	int nrow, ncol;
	int nu;
	double range;
	int cyclic;

public:
	Qmatrix();
	virtual ~Qmatrix();

	int nonzero, N;
	std::vector<std::vector<int> > graph;

	void makeGraph(int *, double *);

	double Q(int , int);
	void lattice2node(int *k, int row, int col);
	void node2lattice(int k, int *row, int *col);
};
int mini(int, int);
int maxi(int, int);
int icmp(const void*, const void*);
//SEXP cMakeQmatrix(SEXP );
#endif /* QMATRIX_H_ */

