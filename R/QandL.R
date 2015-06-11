# Making and solving of the Q matrix, mostly using L
# 
# Author: tarajala
###############################################################################

###############################################################################
## make the Q matrix from scratch using Matern 2D approximation
##
rfhc_makeQ<-function(nrow, ncol, nu, range, cyclic) {
	# using C function
	Qraw<-.External("cMakeQmatrix",
			as.integer(c(ncol, nrow, nu, cyclic)),
			as.numeric(range),
			PACKAGE="rfhc")
	n<-ncol*nrow
	
	# sparse matrix format 
	spMatrix(i=Qraw[[1]]+1, j=Qraw[[2]]+1, x=Qraw[[3]], nrow=n, ncol=n)
}
## full matern cov inversion
#rfhc_makeQ.full<-function(nrow, ncol, nu, range, cyclic){
#	
#	Cf<-function(r, range, nu){
#	kappa <- sqrt(8.0*nu)/range/30
#	(kappa*r)^nu * besselK(kappa*r, nu)/(2^(nu-1)*gamma(nu))	
#	}
#	D<-Cf(extendedgriddists, range, nu)
#	diag(D)<-1.0001
#	forceSymmetric(as(solve(D),"sparseMatrix"))
#}
###############################################################################
## make the Q matrix from scratch using the real deal covariance
## and invert.
##
## rfhc_STATE needs to have the datadists matrix.
rfhc_makeQtrue<-function() {
	Cf<-function(r, range, nu){
		kappa <- sqrt(8.0*nu)/range
		(kappa*r)^nu * besselK(kappa*r, nu)/(2^(nu-1)*gamma(nu))	
	}
	D<-Cf(rfhc_STATE$datadists, rfhc_STATE$range, rfhc_STATE$nu)
	diag(D)<-1
	rfhc_STATE$Qdata<<-solve(D)
}
###############################################################################
## Solving etc with matrices
##
rfhc_SOLVE_LAB<-function(L, B) solve(L, B, system="L")
rfhc_SOLVE_LtAB<-function(L, B) solve(L, B, system="Lt")
rfhc_SOLVE_Lvb<-function(L, b) solve(L, b, system="L")
rfhc_SOLVE_Ltvb<-function(L, b) solve(L, b, system="Lt")
rfhc_SOLVE_QAB<-function(L, B) solve(L, B, system="A")
rfhc_Q2L<-function(Q) Cholesky(Q, LDL=FALSE, super=TRUE, perm=FALSE)
## the perm parameter should be set to TRUE, but need to figure out how the
## UNpermutation then works.
rfhc_SOLVE_INVQ<-function(L)solve(L, system="A")
rfhc_LOGDET<-function(L)determinant(L, log=TRUE)$modulus[1]
###############################################################################
