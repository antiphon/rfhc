#########################################################
##
## run test on the package with function call test()
##
test<-function(...) {
#	simulate_condGMRF(50, 50, 10, 0.2, 2, ...) #ok
	x<-simulate.rfhc(nx=25, r.field.pars=c(-1.5,3,0.2,2), dbg=2) # ok
	fitRFHC(x, ...)
}

simulateGMRF_slow<-function(ncol=50, nrow=50, 
		tau=3, range=0.3, nu=2, expand, 
		xlim=c(0,1), ylim=c(0,1),
		cyclic=FALSE, dbg=0) {
	#
	## if embedding expansion missing expand with 
	## template range:
	if(missing(expand)) expand<-nu*2+1
	#
	## expanded size of the grid
	#
	IO<-rfhc_insideoutside(nrow, ncol, expand)
	#
	## make a state obect which holds important information
	## needed in simulation (to avoid global varialbes)
	#
	rfhc_STATE<<-list(IO=IO, grid=list(nrow=nrow, ncol=ncol, expand=expand), 
			tau=tau, range=range, nu=nu, cyclic=cyclic,
			new=TRUE, save=FALSE, xlim=xlim, ylim=ylim, 
			N=ncol*nrow, Nex=(ncol+2*expand)*(nrow+2*expand), dbg=dbg)
	#
	## simulate using embedding (not true but close?)
	#
	rfhc_sample_GMRF_slow()
	xstep<-diff(xlim)/ncol
	ystep<-diff(ylim)/nrow
	xvec<-seq(xlim[1]+xstep/2, xlim[2]-xstep/2, by=xstep)
	yvec<-seq(ylim[1]+ystep/2, ylim[2]-ystep/2, by=ystep)
	x<-im(matrix(rfhc_STATE$x, byrow=FALSE, nrow), xvec, yvec, xrange=xlim, yrange=ylim)
	x<-rfhc_class_Field(x)
	#
	## save parameters
	#
	x$parameters<-data.frame(tau=tau, range=range, nu=nu,
			cyclic=cyclic,
			nrow=nrow, ncol=ncol, expand=expand)
	x
}

rfhc_sample_GMRF_slow<-function() {
	if(is.null(rfhc_STATE$grid))stop("rfhc_sample_GMRF: Grid not provided in rfhc_STATE, major error.")
	#
	## if IN/OUT node vectors not provided for the extended grid
	if(is.null(rfhc_STATE$IO)) 
		rfhc_STATE$IO<<-rfhc_insideoutside(rfhc_STATE$grid$nrow, rfhc_STATE$grid$ncol, rfhc_STATE$grid$expand)
	#
	## if we want to compute new Q and L
	#
	if(rfhc_STATE$new) {
		#
		## scale the range
		#
		if(is.null(rfhc_STATE$rangescaler))
			rfhc_STATE$rangescaler<<-rfhc_STATE$grid$ncol/diff(rfhc_STATE$xlim)
		#
		## compute Q and L
		#
		ex2<-rfhc_STATE$grid$expand*2
		rfhc_STATE$Q.ex <<- rfhc_makeQ(rfhc_STATE$grid$nrow+ex2, rfhc_STATE$grid$ncol+ex2, rfhc_STATE$nu, rfhc_STATE$range*rfhc_STATE$rangescaler, rfhc_STATE$cyclic)
		rfhc_STATE$L.ex <<- rfhc_Q2L(rfhc_STATE$Q.ex)
		rfhc_STATE$Q    <<- rfhc_STATE$Q.ex[rfhc_STATE$IO$IN, rfhc_STATE$IO$IN]
		rfhc_STATE$L    <<- rfhc_Q2L(rfhc_STATE$Q)
		rfhc_STATE$new  <<- 0
		if(rfhc_STATE$dbg>1)cat("[new Q]")
	}
	#
	z<-rnorm(rfhc_STATE$Nex)
	v<-rfhc_SOLVE_Ltvb(rfhc_STATE$L.ex, z/sqrt( rfhc_STATE$tau) )
	rfhc_STATE$x<<-v[rfhc_STATE$IO$IN]
	#
	## all done.
	#
}
