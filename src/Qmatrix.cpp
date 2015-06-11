/*
 * Qmatrix.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: Tuomas Rajala <tuomas.rajala@jyu.fi>
 */

#include "Qmatrix.h"

Qmatrix::Qmatrix() {
	// TODO Auto-generated constructor stub

}

Qmatrix::~Qmatrix() {
	// TODO Auto-generated destructor stub
}
/**********************************************************
 * Initialise the Matern GMRF 2D lattice neighbourhood
 * system, i.e. the Q matrix. Note that we use tau*Q as precision.
 */

void Qmatrix::makeGraph(int *options, double *theta){
	int nb_row, nb_col, nnb, ir, ic, node1, node2, irow, icol, isnew;
	std::vector<int > neighs;

	ncol=options[0]; //nx
	nrow=options[1]; //ny
	nu=options[2];//g->nu=2;
	cyclic=options[3];//cyclic
	range=theta[0]; //range


	// make the neighbourhood links
	nb_row = 1+nu;
	nb_col = 1+nu;
	nnb = (2*nb_row+1)*(2*nb_col+1);
	N=nrow*ncol;
	graph.resize(N);
	nonzero=0;
	// make links i~j that have Q[i,j]>0 inside square neighbourhood of i.
	if(nnb) {
		for(node1=0; node1<N; node1++){
			neighs.clear();
			node2lattice(node1, &irow, &icol);
			if(cyclic){
				for(ir= irow-nb_row; ir<=irow +nb_row; ir++)
					for(ic = icol -nb_col; ic <= icol+nb_col; ic++){
						lattice2node(&node2, (ir+nrow)%nrow, (ic+ncol)%ncol);
						if(node2!=node1)
//							if(Q(node1, node2)!=0)
								neighs.push_back(node2);
					}
			} else { // not cyclic
				for(ir=maxi(0, irow-nb_row); ir<= mini(nrow-1, irow+nb_row); ir++)
					for(ic = maxi(0, icol-nb_col); ic<=mini(ncol-1, icol+nb_col); ic++){
						lattice2node(&node2, ir, ic);
						if(node2!=node1)
//							if(Q(node1, node2)!=0)
							{
								neighs.push_back(node2);
							}

					}
			}// ok neighbours computed.
			//  Now save to main storage, remove possible duplicates.
			graph.at(node1).clear();
			if(neighs.size()>0){
				graph.at(node1).push_back(neighs.at(0));
				for(ir=1; ir < (int)neighs.size(); ir++) {
					isnew = 1;
					for(ic=0; ic < (int)graph.at(node1).size(); ic++)
						if(neighs.at(ir)==graph.at(node1).at(ic)) {
							isnew=0;
							break;
						}
					if(isnew){
						graph.at(node1).push_back(neighs.at(ir));
					}
				}
			}
			nonzero+= 1 + neighs.size(); // diagonals are also nonzero, even though they are not listed
		}
	} //  Q matrix, as a sparse matrix, computed.
}


/**********************************************************
 * Return Q[i,j] of Matern GMRF 2D lattice neighbourhood matrix.
 * Borrowed from GMRFLib/matern.c, GMRFLib_matern2d(...)
 */
double Qmatrix::Q(int node1, int node2){
	double kappa, a, var=0, val=0;
	kappa=sqrt(8.0*nu)/range;
	a= 4.0+pow(kappa,2.0);

	if(nu>0){
		var=pow(kappa, -2.0*nu)/(4.0*PI*nu);
	}
	else Rprintf("Q: nu=0 not supported, variance set to 0.\n");

	// if Q[i,i]
	if(node1==node2){
		switch (nu) {
		case 0:
			val = a;
			break;
		case 1:
			val = 4.0 + pow(a,2.0);
			break;
		case 2:
			val = a * (pow(a,2.0) + 12.0);
			break;
		case 3:
			val = pow(pow(a,2.0) + 6.0, 2.0) + 12 * pow(a,2.0);
			break;
		default:
			Rprintf("Q: this should not happend.");
		}
		return val * var;
	}
	// i!=j
	int irow, icol, jrow, jcol, drow, dcol, dmin, dmax;

	node2lattice(node1, &irow, &icol);
	node2lattice(node2, &jrow, &jcol);

	drow = abs(irow - jrow);
	dcol = abs(icol - jcol);
	if(cyclic){
		drow = mini(drow, nrow-drow);
		dcol = mini(dcol, ncol-dcol);
	}
	dmax=maxi(drow, dcol);
	dmin=mini(drow, dcol);

	switch(nu){
	case 0:
		if (dmin == 0 && dmax == 1) {
			val = -1.0;
		}
		break;
	case 1:
		switch (dmin) {
		case 0:
			switch (dmax) {
			case 1:
				val = -2.0 * a;
				break;
			case 2:
				val = 1.0;
				break;
			default:
				val = 0.0;
			}
			break;
		case 1:
			val = (dmax == 1 ? 2.0 : 0.0);
		}
		break;
	case 2:
		switch (dmin) {
		case 0:
			switch (dmax) {
			case 1:
				val = -3.0 * (pow(a,2.0) + 3.0);
				break;
			case 2:
				val = 3.0 * a;
				break;
			case 3:
				val = -1.0;
			}
			break;
		case 1:
			switch (dmax) {
			case 1:
				val = 6.0 * a;
				break;
			case 2:
				val = -3.0;
			}
			break;
		}
		break;
	case 3:
		switch (dmin) {
		case 0:
			switch (dmax) {
			case 1:
				val = -4.0 * a * (pow(a,2.0) + 9.0);
				break;
			case 2:
				val = 2.0 * (3.0 * pow(a,2.0) + 8.0);
				break;
			case 3:
				val = -4.0 * a;
				break;
			case 4:
				val = 1.0;
			}
			break;
		case 1:
			switch (dmax) {
			case 1:
				val = 12.0 * (pow(a,2.0) + 2.0);
				break;
			case 2:
				val = -12.0 * a;
				break;
			case 3:
				val = 4.0;
			}
			break;
		case 2:
			val = (dmax == 2 ? 6.0 : 0.0);
			break;
		}
		break;
	default:
		Rprintf("Q: should not happend.\n");
	}
	return val * var;
}

/**********************************************************/

void Qmatrix::lattice2node(int *k, int row, int col){
	*k = row + col*nrow;
}

void Qmatrix::node2lattice(int k, int *row, int *col){
	*col = k/nrow;
	*row = k - (*col)*nrow;
}
/**********************************************************/
int maxi(int x, int y) {
	if(x>y) return x;
	return y;
}
int mini(int x, int y) {
	if(x<y) return x;
	return y;
}

// R call
extern "C"{

SEXP cMakeQmatrix(SEXP Args) {

	int *options;
	double *theta;

	Qmatrix g;
	Args = CDR(Args);
	options = INTEGER(CAR(Args));
	Args = CDR(Args);
	theta = REAL(CAR(Args));
	g.makeGraph(options, theta);

	// Palautetaan SEXP
	SEXP res;
	SEXP ii, jj, v;
	PROTECT(res = allocVector(VECSXP, 3));
	PROTECT(ii = allocVector(INTSXP, g.nonzero));
	PROTECT(jj = allocVector(INTSXP, g.nonzero));
	PROTECT(v = allocVector(REALSXP, g.nonzero));
	int k=0;
	for(int i=0;i<g.N; i++){
		INTEGER(ii)[k] = i;
		INTEGER(jj)[k] = i;
		REAL(v)[k]= g.Q(i,i);
		k++;
		for(int j=0; j < (int)g.graph.at(i).size(); j++) {
			INTEGER(ii)[k] = i;
			INTEGER(jj)[k] = g.graph.at(i).at(j);
			REAL(v)[k]= g.Q(i, g.graph.at(i).at(j));
			k++;
		}
	}
	SET_VECTOR_ELT(res, 0, ii);
	SET_VECTOR_ELT(res, 1, jj);
	SET_VECTOR_ELT(res, 2, v);
	UNPROTECT(4);

	return res;
}

} // extern

