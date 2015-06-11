/*
 * Grid.cpp
 *
 *  Created on: 7.7.2009
 *      Author: tarajala
 */

#include "Grid.h"


Grid::Grid() {
	// TODO Auto-generated constructor stub

}

Grid::~Grid() {
	// TODO Auto-generated destructor stub
}
void Grid::Init(SEXP X)// from image-object
{
	int i,j;
	double *v;
	this->xlim = REAL(getListElement(X, "xrange"));
	this->ylim = REAL(getListElement(X, "yrange"));
	this->ny = INTEGER(getListElement(X, "dim"))[0];
	this->nx = INTEGER(getListElement(X, "dim"))[1];
	this->xstep = REAL(getListElement(X,"xstep"))[0];
	this->ystep = REAL(getListElement(X,"ystep"))[0];
	this->n = nx*ny;
	this->dxdy = xstep*ystep;
	this->values.resize(this->n, 0.0);
	this->values2.resize(this->n, 0.0);
	this->bincount.resize(this->n, 0);
	v = REAL(getListElement(X, "v"));
	integral =0;
	for(i=0; i<this->ny;i++){
		for(j=0; j<nx;j++){
			values.at(i*nx+j) = v[i+j*ny];
			integral = integral + values.at(i*nx+j)*dxdy;
		}
	}
	this->getBoxWeightsp = &Grid::getBoxWeightsAll1;

}

SEXP Grid::toSEXP()
{
	SEXP res;
	SEXP *vecz;
	double *pz;
	int i;
	PROTECT(res = allocVector(VECSXP, 1));
	vecz = new SEXP;
	PROTECT(*vecz = allocVector(REALSXP, this->n));
	pz =    REAL(*vecz);
	for(i=0;i<this->n;i++)
	{
		pz[i] = this->getValue(&i);
	}
	SET_VECTOR_ELT(res, 0, *vecz);
	UNPROTECT(2);
	return res;

}

int Grid::getSize(int *i)
{
	if(*i == 0 ) return this->nx;
	else if(*i == 1 ) return this->ny;
	return this->n;
}

int Grid::getSize()
{
	return this->n;
}

void Grid::setValue(int *i, int *j, double *v)// i = row, j = col, top-bottom, left-right
{
	if((*j >= nx) | (*j <0) | (*i<0) | (*i >= ny));
	else{
		this->values.at(  (*i)*(this->nx)+(*j) ) = *v;
	}
}

void Grid::setValue(int *i, double *v)// i = row, j = col, top-bottom, left-right
{
	if((*i<0) | (*i >= n));
	else{
		this->values.at( *i ) = *v;
	}
}

double Grid::getValue(int *i, int *j)// i = row, j = col, top-bottom, left-right
{
	if((*j >= nx) | (*j <0) | (*i<0) | (*i >= ny)) return -999999;
	return this->values.at(  (*i)*(this->nx)+(*j) );
}

double Grid::getValue(int *i)// i = number in the vector
{
	if((*i<0) | (*i >= this->n)) return -999999;
	return this->values.at(*i);
}

double Grid::getValue(double x, double y)// (x,y) real location, get nearest grid value
{
	if((x<xlim[0]) | (x >=xlim[1]) | (y<ylim[0]) | (y>= ylim[1])) {
		Rprintf("[Grid: requested value outside window (x=%f, y=%f)]",x, y);
		return 0;
	}
	int i,j;
	i = (int) floor((y-ylim[0])/this->ystep);
	j = (int) floor((x-xlim[0])/this->xstep);
	return this->getValue(&i, &j);
}

void Grid::addValue(int *i, int *j, double *v)// i = row, j = col, top-bottom, left-right
{

	if((*j < nx) & (*j >= 0) & (*i >= 0) & (*i < ny))
	{
		values[ (*i)*(nx)+(*j) ] = values.at( (*i)*(nx)+(*j) ) + (*v);
		bincount[(*i)*(nx)+(*j)]++;
	}
}
void Grid::multiplyValue(int *i, int *j, double *v)// i = row, j = col, top-bottom, left-right
{

	if((*j < nx) & (*j >= 0) & (*i >= 0) & (*i < ny))
	{
		values[ (*i)*(nx)+(*j) ] = values.at( (*i)*(nx)+(*j) ) * (*v);
		bincount[(*i)*(nx)+(*j)]++;
	}
}



void Grid::addValue2(int *i, int *j, double *v)// i = row, j = col, top-bottom, left-right
{
	if((*j < nx) & (*j >= 0) & (*i >= 0) & (*i < ny))
	{
		values2[ (*i)*(nx)+(*j) ] = values2.at( (*i)*(nx)+(*j) ) + (*v);
	}
}

void Grid::addValueXY(double *x, double *y, double *v)// x,y coords, aggregate
{
	int i,j;
	i = (int) floor((*y-ylim[0])/this->ystep);
	j = (int) floor((*x-xlim[0])/this->xstep);
	this->addValue(&i, &j, v);
}
void Grid::addValue2XY(double *x, double *y, double *v)// x,y coords, aggregate
{
	int i,j;
	i = (int) floor((*y-ylim[0])/this->ystep);
	j = (int) floor((*x-xlim[0])/this->xstep);
	this->addValue2(&i, &j, v);
}

void Grid::getCenter(int *i, int *j, double *c)
{
	c[0] = ((double)(*j)+0.5) * this->xstep + this->xlim[0];
	c[1] = ((double)(*i)+0.5) * this->ystep + this->ylim[0];
}

void Grid::setIntegral(double *v){
	integral = *v;
}
void Grid::addIntegral(double *v){
	integral = integral + dxdy*(*v);
}

double Grid::getIntegral(){
	return integral;
}

double Grid::getXstep()
{
	return this->xstep;
}

double Grid::getYstep()
{
	return this->ystep;
}

void Grid::takeMean()
{
	int i;
	for(i=0;i < this->n;i++)
		if(this->values2.at(i)>0)
			this->values.at(i) = this->values.at(i)/this->values2.at(i);
}

void Grid::calcBoxWeights(double h)
{
//  edge correction weights for smoothing.
//	first version, not optimised. Direct geometrical computation.
	int i,j;
	this->weights.clear();
	double dx, dy, d, c[2], a1, a2, b1, b2, u, tx1, tx2, ty1, ty2, wij;

	for(i=0; i<ny; i++)
		for(j=0; j<nx; j++){
			this->getCenter(&i, &j, c);
			dx=fmin(h, fmin(c[0]-xlim[0], xlim[1]-c[0] ));
			dy=fmin(h, fmin(c[1]-ylim[0], ylim[1]-c[1] ));
			if(fmin(dx, dy)<h){ // wheter the circle touches the border
				a1=a2=acos(dx/h);
				b1=b2=acos(dy/h);
				u=sqrt(dx*dx+dy*dy);
				if(u<h){
					a2=acos(dx/u);
					b2=acos(dy/u);
				}
				d=2*PI-a1-a2-b1-b2;
				tx1=tan(a1)*dx;
				tx2=tan(a2)*dx;
				ty1=tan(b1)*dy;
				ty2=tan(b2)*dy;
				wij=1.0/(((d/2.0)*h*h+dx*(tx1+tx2)/2.0+dy*(ty1+ty2)/2.0)/(h*h*PI));
				this->weights.push_back(wij);
			}
			else{
				this->weights.push_back(1.0);
			}
		}


	this->getBoxWeightsp = &Grid::getBoxWeightsComputed;
}

double Grid::getBoxWeights(int *i, int *j)
{
	return (this->*getBoxWeightsp)(i, j);
}

double Grid::getBoxWeightsAll1(int *i, int *j)
{
	return 1.0;
}

double Grid::getBoxWeightsComputed(int *i, int *j)
{
	return weights.at((*i)*nx+(*j));
}

void Grid::smoothen()
{
//	int i,j;
//	for(i=0;i < this->ny-1;i++) // first row
	Rprintf("smoother not yet implemented");
}



/*******************************************************************/
std::vector<int> Grid::getNeighbours(int *i)
{
	int r,c;
	std::vector<int> x;
	r = (int) (double)*i/(double)nx;
	c = *i - r*nx;
	x.clear();
	if(r == 0)
	{
		if(c==0)
		{
			x.push_back(1);
			x.push_back(nx);
			x.push_back(nx+1);
		}
		else
		{
			if(c==(nx-1) )
			{
				x.push_back(nx-2);
				x.push_back(nx+nx-1);
				x.push_back(nx+nx-2);
			}
			else
			{
				x.push_back(*i-1);x.push_back(*i+1);
				x.push_back(*i-1 + nx);x.push_back(*i + nx);x.push_back(*i+1 + nx);

			}
		}
	}
	else if(r == (ny-1) )
	{
		if(c==0)
		{
			x.push_back(n-nx+1);
			x.push_back(n-nx-nx);
			x.push_back(n-nx-nx+1);
		}
		else
		{
			if(c==(nx-1) )
			{
				x.push_back(n-2);
				x.push_back(n-nx-1);
				x.push_back(n-nx-2);
			}
			else
			{
				x.push_back((n-nx)+c-1);x.push_back((n-nx)+c+1);
				x.push_back((n-nx-nx)+c-1);x.push_back((n-nx-nx)+c);x.push_back((n-nx-nx)+c+1);
			}
		}
	}
	else if(c == 0)
	{
		x.push_back((r-1)*nx);x.push_back((r-1)*nx+1);
		x.push_back(r*nx+1);
		x.push_back((r+1)*nx);x.push_back((r+1)*nx+1);
	}
	else if(c == (nx-1) )
	{
		x.push_back(r*nx-1);x.push_back(r*nx-2);
		x.push_back((r+1)*nx-2);
		x.push_back((r+2)*nx-2);x.push_back((r+2)*nx-1);
	}
	else
	{
		x.push_back((r-1)*nx+c-1);x.push_back((r-1)*nx+c);x.push_back((r-1)*nx+c+1);
		x.push_back(r*nx+c-1);x.push_back(r*nx+c+1);
		x.push_back((r+1)*nx+c-1);x.push_back((r+1)*nx+c);x.push_back((r+1)*nx+c+1);
	}

	return x;
}

