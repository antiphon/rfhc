/*
 * SimRFHC.cpp
 *
 *  Created on: Oct 12, 2011
 *      Author: antiphon
 */

#include "SimRFHC.h"

SimRFHC::SimRFHC() {
	// TODO Auto-generated constructor stub

}

SimRFHC::~SimRFHC() {
	// TODO Auto-generated destructor stub
}


void SimRFHC::Init(SEXP pp0, SEXP rf, int *d0, double *pars){
	pp.Init(pp0);
	debug = *d0;
	//	pars: c(as.numeric(n), giveup, gamma, iter),
	N      = (int)pars[0];
	giveup = (int)pars[1];
	gamma  =      0;
	iter   = (int)pars[2];
	typeOfInteraction = (int) pars[3];
	numberAIest = 0;//(int) pars[5];

	if(typeOfInteraction==1){
		checkInteractionp = &SimRFHC::checkHC;
		logdensityp = &SimRFHC::logdensityHC;
	}
	else{
		checkInteractionp = &SimRFHC::checkHS;
		logdensityp = &SimRFHC::logdensityHS;
	}
	if(debug>1)Rprintf("Got[N=%i, giveup=%i, iter=%i, type=%i]\n", N,giveup,iter,typeOfInteraction);
	rfield.Init(rf);
	double m;
	for(int i=0; i<pp.size(); i++) {
		m = rfield.getValue(pp.getX(i), pp.getY(i));
		pp.setMass(&i, &m);
	}
	paphist = new double[iter];
}
/********************************************************/
/********************************************************/
void SimRFHC::Run(){
//	generate first a starting pattern:
	int t=0, i, j, fails=0, acc;
	bool newok;
	double *xy, m, logden;
//	double p, logpap0, logpap1;
	xy = new double[2];
	int numbers[2]; numbers[0]=0; numbers[1]=0;
	/*
	  1. generate an acceptable starting pattern
	 */
	j=-1;
	while((pp.size() < N) & (t < giveup)) {
		runifpoint(xy);
		m = rfield.getValue(xy[0], xy[1]);
		// check hard spheres
		newok = checkInteraction(xy, &m, &j);
		if(newok){
			// no violation, add to pattern
			pp.addPoint(xy[0], xy[1]);
			i = pp.size()-1;
			pp.setMass(&i, &m);
			numbers[0]++;
			fails=0;
		}
		else{
			fails++;
			if((fails > 1000) & (pp.size()>0)){
				i = (int) ((double)pp.size())*unif_rand();
				pp.deletePoint(i);
				numbers[1]++;
				fails = 0;
			}
		}
		t++;
		if(debug)Rprintf("\r[n=%i, b=%i, d=%i]    ", pp.size(), numbers[0], numbers[1]);
	}
	if((giveup <= t) & (debug>0)) Rprintf("Gave up after %i attempts.\n", giveup);

	// compute the non-normalised density value:
	logden=0;
	logden = logdensity();
	/*
		2. shuffle using MH algorithm:
	*/
	t=0;
	acc=0;
	while(t < iter) {
		runifpoint(xy);
		m = rfield.getValue(xy[0], xy[1]);
		j = (int) ((double)pp.size())*unif_rand();
		// check violations
		newok = checkInteraction(xy, &m, &j);


		if(newok){ // no hard sphere violation, flip coin:
//			logpap0 = -gamma*PI*pow(pp.getMass(&j),2);
//			logpap1 = -gamma*PI*pow(m,2);
//			p = fmin(exp(logpap1 -logpap0),1);
//			 if we update:
//			if(unif_rand() < p ){
//				logden = logden - logpap0 + logpap1;
				pp.deletePoint(j);
				pp.addPoint(xy[0], xy[1]);
				i = pp.size()-1;
				pp.setMass(&i, &m);
				acc++;
//			}
		}
//		paphist[t] = logden;
		t++;
		if(debug)Rprintf("\r[%i/%i, acc=%i]", t, iter,acc);
	}
	delete xy;
}

/********************************************************/
void SimRFHC::runifpoint(double *xy) {
	xy[0] = unif_rand() * (pp.xlim[1]-pp.xlim[0]) + pp.xlim[0];
	xy[1] = unif_rand() * (pp.ylim[1]-pp.ylim[0]) + pp.ylim[0];
}
/*******************************************************/


/*******************************************************
	The interactions: HC and HS
*******************************************************/

bool SimRFHC::checkInteraction(double *xy, double *m, int *j){
	return (this->*checkInteractionp)(xy, m, j);
}
/********************************************************/
bool SimRFHC::checkHC(double *xy, double *m, int *j){
	double mi, d;
	for(int i=0; i<pp.size(); i++) {
		if(i!=(*j)){
			d = sqrt(pow(xy[0]-pp.getX(i),2.0)+pow(xy[1]-pp.getY(i),2.0));
			mi = pp.getMass(&i);
			if(d < fmax(*m,mi) ) return false;
		}
	}
	return true;
}
/********************************************************/
bool SimRFHC::checkHS(double *xy, double *m, int *j){
	double mi, d;
	for(int i=0; i<pp.size(); i++) {
		if(i!=(*j)){
			d = sqrt(pow(xy[0]-pp.getX(i),2.0)+pow(xy[1]-pp.getY(i),2.0));
			mi = pp.getMass(&i);
			if(d < (*m+mi) ) return false;
		}
	}
	return true;
}
/********************************************************/
/********************************************************/

/********************************************************
    The log densities if we want to incorporate
    AI component: full pattern
********************************************************/
double SimRFHC::logdensity() {
	return (this->*logdensityp)();
}
/********************************************************
 * Log-density of Hard Core model:
 * need to compute the total area covered by the
 * discs b(xi, Ri).
 *
 * This version estimates the area with
 * uniformly random points.
 *
 */
double SimRFHC::logdensityHC(){
	double xy[2];
	int i, j, hit=0;
	for(i=0; i<numberAIest; i++){
		runifpoint(xy);
		for(j=0;j<pp.size();j++){
			if((pow(xy[0]-pp.getX(j),2)+pow(xy[1]-pp.getY(j),2))<pow(pp.getMass(&j),2)){
				hit++;
				break;
			}
		}
	}
	// compute the covered area and multiply with -gamma
	return -gamma*((double)hit/(double)numberAIest)*pp.getWindowArea();
}
/****************************************************
 * Log-density of Hard Spheres model:
 * this is simply the sum of areas as the discs
 * b(xi, Ri) dont overlap.
 */
double SimRFHC::logdensityHS(){
	double logden=0;
	for(int i=0; i<pp.size();i++)
		logden = logden - pow(pp.getMass(&i),2);
	return gamma*PI*logden;
}

/********************************************************/
/********************************************************/
/********************************************************/
SEXP SimRFHC::toSEXP(){
	SEXP res, *x, *y, *m, *hist;
	double *px, *py, *pm, *phist;
	int i;
	x = new SEXP;
	y = new SEXP;
	m = new SEXP;
	PROTECT(*x = allocVector(REALSXP, pp.size()));
	PROTECT(*y = allocVector(REALSXP, pp.size()));
	PROTECT(*m = allocVector(REALSXP, pp.size()));
	px = REAL(*x);
	py = REAL(*y);
	pm = REAL(*m);
	for(i=0; i<pp.size();i++){
		px[i]=pp.getX(i);
		py[i]=pp.getY(i);
		pm[i]=pp.getMass(&i);
	}

	hist = new SEXP;
	PROTECT(*hist = allocVector(REALSXP, iter));
	phist=REAL(*hist);
//	memcpy(phist, paphist, iter);
	for(i=0;i<iter;i++)phist[i]=paphist[i];

	PROTECT(res = allocVector(VECSXP, 4));
	SET_VECTOR_ELT(res, 0, *x);
	SET_VECTOR_ELT(res, 1, *y);
	SET_VECTOR_ELT(res, 2, *m);
	SET_VECTOR_ELT(res, 3, *hist);
	UNPROTECT(5);
	delete x;
	delete y;
	delete m;
	delete hist;
	delete paphist;
	return res;
}
/********************************************************/
/********************************************************/
/********************************************************/
extern "C" {

SEXP simRFHC_c(SEXP Args) {
	SEXP pp, rfield;
	SimRFHC simulator;
	int *dbg;
	double *pars;

	Args = CDR(Args);
	dbg = INTEGER(CAR(Args)); // if debug messages

	Args = CDR(Args);
	pp = CAR(Args); // init pp

	Args = CDR(Args);
	rfield = CAR(Args); // init r-field

	Args = CDR(Args);
	pars = REAL(CAR(Args)); // parameters, only n a.t.m.

//	init
	simulator.Init(pp, rfield, dbg, pars);

	GetRNGstate();
	simulator.Run();
	PutRNGstate();
	return simulator.toSEXP();
}



}
