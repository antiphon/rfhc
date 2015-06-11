#include "Pp.h"
/********************************************************************************************/
Pp::Pp()
{
}
/********************************************************************************************/
Pp::~Pp()
{
}
/********************************************************************************************/
void Pp::Init(SEXP Argspp)
{
	int i, *type;
	double *x, *y, *m, *m2;
	Point *p;
	tor = 0;
	toroidal = &tor;
	distances_calculated=0;
	this->distancep = &Pp::distance_calc;
	this->nndistancep = &Pp::nndistance_calc;
	this->m = length(getListElement(Argspp, "x"));
	x = REAL(getListElement(Argspp, "x"));
	y = REAL(getListElement(Argspp, "y"));
	type = INTEGER(getListElement(Argspp, "types"));
	m = REAL(getListElement(Argspp, "mass"));
	m2 = REAL(getListElement(Argspp, "mass2"));

	windowArea =  REAL(getListElement(Argspp, "area"))[0];


	points.clear();
	for(i=0; i < this->size(); i++)
	{
		p = new Point(x[i],y[i],type[i]);
		p->setId(&i);
		p->setMass(&(m[i]));
		p->setMass2(&(m2[i]));
		points.push_back(*p);

	}
	this->m = points.size();
	this->ntypes = type[this->size()];

	xlim = REAL(getListElement(getListElement(Argspp, "window") ,"x"));
	ylim = REAL(getListElement(getListElement(Argspp, "window") ,"y"));
	zlim = new double [2];zlim[0]=0; zlim[1]=0;
}

/********************************************************************************************/
double  Pp::getX(int i) {return this->points[i].getX();}
double  Pp::getY(int i) {return this->points[i].getY();}
double  Pp::getZ(int i) {return this->points[i].getZ();}
int     Pp::getT(int i) {return this->points[i].getT();}
int 	Pp::size()      {return this->m;   }
int     Pp::nsize(int *i){return this->points[*i].nsize();}
int     Pp::getCluster(int i){return this->points[i].getCluster();}
double  Pp::getMass(int *i){return this->points[*i].getMass();}
void    Pp::setMass(int *i, double *x){this->points[*i].setMass(x);}
double  Pp::getMass2(int *i){return this->points[*i].getMass2();}
void    Pp::setMass2(int *i, double *x){this->points[*i].setMass2(x);}
int     Pp::getNtypes(){return this->ntypes;}
double Pp::getWindowArea(){return this->windowArea;}
/********************************************************************************************/
void    Pp::movePoint(int *i, double *x, double *y){this->points[*i].move(x, y);}
void    Pp::clearNeighbourhood(int *i) {this->points[*i].clearNeighbourhood();}
/********************************************************************************************/

void Pp::distance_precalc()
{
	int i,j;
	double d;
	distMatrix.resize(m);
	for(i=0;i<m;i++)distMatrix.at(i).resize(m);
	for(i=0; i<m-1;i++){
		for(j=i;j<m;j++){
			d=distance_calc(&i, &j);
			distMatrix.at(i).at(j)=d;
			distMatrix.at(j).at(i)=d;
		}
	}
	distances_calculated = 1;
	this->distancep = &Pp::distance_fetch;
}

/********************************************************************************************/
double Pp::distance(int *i, int *j)
{
	return (this->*distancep)(i, j);
}
/********************************************************************************************/
double Pp::distance_fetch(int *i, int *j)
{
	return distMatrix.at(*i).at(*j);
}
/********************************************************************************************/
double Pp::distance_calc(int *i, int *j)
{
	if(*i==*j) return 0.0;
	if(*this->toroidal)
		return	sqrt(
					pow( fminf( this->xlim[1]-this->xlim[0]-fabs(getX(*i)-getX(*j)) , fabs(getX(*i)-getX(*j)) ) ,2.0) +
					pow( fminf( this->ylim[1]-this->ylim[0]-fabs(getY(*i)-getY(*j)) , fabs(getY(*i)-getY(*j)) ) ,2.0) +
					pow( fminf( this->zlim[1]-this->zlim[0]-fabs(getZ(*i)-getZ(*j)) , fabs(getZ(*i)-getZ(*j)) ) ,2.0)   );
	else
		return 	sqrt(
				pow( getX(*i)- getX(*j)  ,2.0) +
				pow( getY(*i)- getY(*j)  ,2.0) +
				pow( getZ(*i)- getZ(*j)  ,2.0)   );

}
/********************************************************************************************/
void Pp::distance_update(int *i)
{
	if(distances_calculated==0)distance_precalc();
	double d;
	int j;
	for(j=0; j<m;j++) {
		d = distance_calc(i, &j);
		distMatrix.at(*i).at(j)=d;
		distMatrix.at(j).at(*i)=d;
	}
}

/********************************************************************************************/
/********************************************************************************************/
void Pp::nndistance_precalc()
{
	int i;
	double nnd;
	nndistvector.resize(m);
	for(i=0; i<m;i++){
		nnd=nndistance_calc(&i);
		nndistvector.at(i)=nnd;
	}
	nndistances_calculated = 1;
	this->nndistancep = &Pp::nndistance_fetch;
}
/********************************************************************************************/
double Pp::nndistance(int *i)
{
	return (this->*nndistancep)(i);
}
/********************************************************************************************/
double Pp::nndistance_calc(int *i)
{
	double *dists = new double[m];
	double d;
	for(int j=0;j<m;j++)dists[j]=distance(i, &j);
	qsort( dists, m, sizeof(double), compare_doubles);
	d=dists[1];
	delete [] dists;
	return d;

}
/********************************************************************************************/
double Pp::nndistance_fetch(int *i)
{
	return nndistvector.at(*i);
}
/********************************************************************************************/
/********************************************************************************************/
void Pp::addPoint(double x, double y)
{
	int i = points.size()+1,t=0;
	double z=0;
	Point *p;
	p = new Point(x, y, t);
	p->setId(&i);
	p->setMass(&z);
	p->setMass2(&z);
	points.push_back(*p);
	this->m++;
	delete p;
}
/********************************************************************************************/
void Pp::deletePoint(int i){
	points.erase(points.begin()+i);
	m--;
}
/********************************************************************************************/
void Pp::addNeighbour(int *i, int *j)
{
	this->points.at(*i).addNeighbour(j);
}
/********************************************************************************************/
void Pp::removeNeighbour(int *i, int *j)
{
	this->points.at(*i).getNeighbour(j);
}
/********************************************************************************************/
int Pp::getNeighbour(int *i, int *j)
{
	return this->points.at(*i).getNeighbour(j);
}
/********************************************************************************************/
void Pp::setCluster(int *i, int *j)
{
	this->points.at(*i).setCluster(j);
}

/**********************************************************************************/
int compare_doubles(const void *a, const void *b)
{
  const double *da = (const double *) a;
  const double *db = (const double *) b;

  return (*da > *db) - (*da < *db);
}

// EOF
