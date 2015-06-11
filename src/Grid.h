/*
 * Grid.h
 *
 *  Created on: 7.7.2009
 *      Author: tarajala
 */

#include <R.h>
#include <Rinternals.h>
#include "Rextras.h"
#include <vector>
#ifndef GRID_H_
#define GRID_H_

class Grid {
	double xstep;
	double ystep;
	double *xlim;
	double *ylim;
	double integral;
	double integral_chemical;
	std::vector<double> values;
	std::vector<double> values2;
	std::vector<int> bincount;
	std::vector<double> weights;
	double (Grid::*getBoxWeightsp)(int*, int*);
	double getBoxWeightsAll1(int*, int*);
	double getBoxWeightsComputed(int*, int*);

public:
	int nx;
	int ny;
	int n;
	double dxdy;

	Grid();
	virtual ~Grid();
	void Init(SEXP);
	int getSize(int *);
	int getSize();
	void   setValue(int *, int *, double *);
	void   setValue(int *, double *);
	double getValue(int *, int *);
	double getValue(int *);
	double getValue(double x, double y);
	void   addValue(int *, int *, double *);
	void   multiplyValue(int *, int *, double *);
	void   addValue2(int *, int *, double *);
	void   addValueXY(double *, double *, double *);
	void   addValue2XY(double *, double *, double *);
	void   getCenter(int *, int *, double *);
	void	setIntegral(double *);
	void   addIntegral(double *);
	double getIntegral();
	double getXstep();
	double getYstep();
	void   takeMean();
	void   smoothen();
	void   calcBoxWeights(double);
	double getBoxWeights(int*, int*);
	int getNsize(int *);
	std::vector<int> getNeighbours(int *);
	SEXP toSEXP();
};

#endif /* GRID_H_ */
