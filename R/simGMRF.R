# Functions for simulating GMRF, conditional and unconditional.
# 
# Author: tarajala
###############################################################################

###############################################################################
### FUNCTIONS TO BE USED DIRECTLY
###

#################################
## simulate unconditional. mean 0
simulateGMRF<-function(ncol=50, nrow=50, 
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
	rfhc_sample_GMRF()
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
#################################
## simulate conditional. mean 0
simulate_condGMRF<-function(ncol=50, nrow=50, 
		tau=5, range=0.2, nu=2, expand, 
		xlim=c(0,1), ylim=c(0,1),
		data=list(x=.5, y=.5, v=0),
		cyclic=FALSE, dbg=0, use.RandomFields=FALSE, ...) {
	#
	## if embedding expansion missing expand with 
	## template range:
	if(missing(expand)) expand<-nu*2+1
	#
	## expanded size of the grid
	#
	IO<-rfhc_insideoutside(nrow, ncol, expand)
	#
	## make the conditional locations to grid nodes
	#
	xcol<-seq(xlim[1], xlim[2], length=ncol)
	yrow<-seq(ylim[1], ylim[2], length=nrow)
	fixed<-rfhc_xy2node(data$x, data$y, xcol, yrow)
	## make a state obect which holds important information
	## needed in simulation 
	#
	rfhc_STATE<<-list(IO=IO, grid=list(nrow=nrow, ncol=ncol, expand=expand), 
			tau=tau, range=range, nu=nu, cyclic=cyclic,
			new=TRUE, save=FALSE, xlim=xlim, ylim=ylim, 
			N=ncol*nrow, Nex=(ncol+2*expand)*(nrow+2*expand),
			x.cond=data$v, fixed=fixed,
			cond.is.new=TRUE, dbg=dbg
			)
	#
	## simulate using embedding (not true but close?)
	#
	## consider also using RandomFields
    if(use.RandomFields){
        rfhc_STATE$data<<-list(x=data$x, y=data$y)
        rfhc_sample_conditional_GRF(...)
    }
    else rfhc_sample_conditional_GMRF()
	
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
#	list(x=x, s=rfhc_STATE)
	x
}

#####################################################################
#####################################################################
### INTERNAL FUNCTIONS
##
#####################################################################
## Main work horse for sampling from Conditional 
## normal random field using Matern 2D 
## Markov approximation.
## 
## MEAN == 0
##
## Uses internal representation of the problem as parameters.
##
rfhc_sample_conditional_GMRF<-function() {
	#
	## Note that conditioning should be in nodeworld,
	# and rfhc_STATE$fixed[i] if sample$x.cond[i] is 
	# to be fixed.
	#
	if(is.null(rfhc_STATE$fixed))stop("rfhc_STATE$fixed missing")
	##
	## start by simulating the unconditional sample:
	#
	rfhc_sample_GMRF()
 	
	#
	## then correct for conditioning:
	#
	if(rfhc_STATE$cond.is.new){
		nf<-length(rfhc_STATE$fixed)
		rfhc_STATE$A   <<- spMatrix(j=rfhc_STATE$fixed, i=1:nf, x=rep(1, nf), nf, rfhc_STATE$N)
		rfhc_STATE$V   <<- rfhc_SOLVE_QAB( rfhc_STATE$L, t(rfhc_STATE$A)/sqrt(rfhc_STATE$tau) )
		rfhc_STATE$W   <<- rfhc_STATE$A%*%rfhc_STATE$V
		rfhc_STATE$Winv<<- solve(rfhc_STATE$W)
		rfhc_STATE$U   <<- rfhc_STATE$Winv%*%t(rfhc_STATE$V)
		rfhc_STATE$Ax  <<- rfhc_STATE$A%*%rfhc_STATE$x
		rfhc_STATE$co  <<- rfhc_STATE$Ax-rfhc_STATE$x.cond

	}
	rfhc_STATE$x<<-as.numeric( rfhc_STATE$x - t(rfhc_STATE$U)%*%rfhc_STATE$co ) 
	# ok
	## for repeated samples
	#
	rfhc_STATE$cond.is.new<<-FALSE
	#
}

#####################################################################
#####################################################################
## Main work horse for sampling from UNCOnditional 
## normal random field using Matern 2D 
## Markov approximation.
##
## Uses internal representation of the problem as parameters.
rfhc_sample_GMRF<-function() {
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
	## simulate: note that Q does not have tau in it
	#
	z<-rnorm(rfhc_STATE$Nex)
	v<-rfhc_SOLVE_Ltvb(rfhc_STATE$L.ex, z/sqrt(rfhc_STATE$tau) )
	rfhc_STATE$x<<-v[rfhc_STATE$IO$IN]
	#
	## all done.
	#
}
#####################################################################
#####################################################################
### do the conditional simulation using RandomFields:
rfhc_sample_conditional_GRF <- function(...){
    require(RandomFields)
    ## build required vectors
    xcol <- seq(rfhc_STATE$xlim[1], rfhc_STATE$xlim[2], length=rfhc_STATE$grid$ncol)
    yrow <- seq(rfhc_STATE$ylim[1], rfhc_STATE$ylim[2], length=rfhc_STATE$grid$nrow)
    dataxy <- cbind(rfhc_STATE$data$x, rfhc_STATE$data$y)
    datav  <- rfhc_STATE$x.cond
    nu     <- rfhc_STATE$nu
    ## match parameters:
    s2    <- 1/rfhc_STATE$tau
    kappa <- sqrt(8.0*nu)/rfhc_STATE$range
    ## simple kriging as we set mu==0
    v<-CondSimu("S", grid=TRUE, x=xcol, y=yrow, given=dataxy, data=datav, 
                model="matern", param=c(0, s2, 0, 1/kappa, nu), pch="", ...)
    rfhc_STATE$x<<-c(t(v))
}
#####################################################################
#####################################################################
#eof






