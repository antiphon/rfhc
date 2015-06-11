#include <R.h>
#include <Rinternals.h>
#include "Point.h"
#include "Rextras.h"
#ifndef PP_H_
#define PP_H_

class Pp
{
	std::vector<Point> points;
	int m;
	int ntypes;
	int tor;
	int *toroidal;
	int distances_calculated, nndistances_calculated;
	double windowArea;
	std::vector<std::vector<double> >distMatrix;
	double (Pp::*distancep)(int*, int*);
	double (Pp::*nndistancep)(int*);
	std::vector<double> nndistvector;

public:
	double *xlim;
	double *ylim;
	double *zlim;


	Pp();
	virtual ~Pp();

	void Init(SEXP);

	double getX(int);
	double getY(int);
	double getZ(int);
	int    getT(int);
	int    getCluster(int);
	int    size();
	int    nsize(int*);
	void   setMass(int *, double *);
	double getMass(int *);
	void   setMass2(int *, double *);
	double getMass2(int *);
	int    getNtypes();
	double getWindowArea();
	void   distance_precalc();
	double distance(int *, int *);
	double distance_calc(int *, int *);
	double distance_fetch(int *, int *);
	void   distance_update(int *);
	void   nndistance_precalc();
	double nndistance(int *);
	double nndistance_calc(int *);
	double nndistance_fetch(int *);
	void addPoint(double, double);
	void deletePoint(int);
	void addNeighbour(int *i, int *j);
	void removeNeighbour(int *i, int *j);
	void clearNeighbourhood(int *i);
	int  getNeighbour(int *i, int *j);
	void setCluster(int *i, int *j);
	void movePoint(int *, double *, double *);
};

#endif /*PP_H_*/
int compare_doubles(const void *a, const void *b);
