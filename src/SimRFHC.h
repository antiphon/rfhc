/*
 * SimRFHC.h
 *
 *  Created on: Oct 12, 2011
 *      Author: Tuomas Rajala <tuomas.rajala@jyu.fi>
 */

#ifndef SIMRFHC_H_
#define SIMRFHC_H_
#include <R.h>
#include <Rinternals.h>
#include "Pp.h"
#include "Grid.h"
#include <vector>
class SimRFHC {
	Pp pp;
	Grid rfield;
	int debug, N, iter, giveup, typeOfInteraction, numberAIest;
	double *paphist, gamma;

	/********************************************************/
	bool (SimRFHC::*checkInteractionp)(double*, double*, int*);
	bool checkInteraction(double *, double *, int *);
	bool checkHS(double *, double *, int *);
	bool checkHC(double *, double *, int *);
	/********************************************************/
	double (SimRFHC::*logdensityp)();
	double logdensity();
	double logdensityHC();
	double logdensityHS();
	/********************************************************/
public:
	SimRFHC();
	virtual ~SimRFHC();
	void Init(SEXP, SEXP, int*, double*);
	void runifpoint(double *);
	void Run();
	SEXP toSEXP();
};

#endif /* SIMRFHC_H_ */

//SEXP simRFHC(SEXP Args);
