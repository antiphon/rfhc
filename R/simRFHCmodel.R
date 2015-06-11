# Direct functions for simulating the varying range hard-core model
# 
# Author: antiphon
###############################################################################

simulateRFHC<-function(n=100,
		iter=1e5,
		ncol=50,
		window=square(),
		r.field.pars=c(mean=-2, tau=5, range=0.4, nu=2),
		start.pattern=NULL,
		pre.r.field=NULL,
		giveup=1e6,
		type="HC",
		dbg=0, seed, cyclic, expand) {
	
	T0<-Sys.time()
	#
	## keep aspect ratio 1 in grid
	#
	nrow<-round(ncol*diff(window$yrange)/diff(window$xrange))
	#
	## seed
	#
	if(missing(seed)){runif(1); seed<-as.integer(format(Sys.time(),"%s"))*runif(1)}
	set.seed(seed)
	#
	## starting pattern
	#
	if(is.null(start.pattern))xy<-rpoint(n=1, win=window)
	else xy<-ppp(start.pattern$x, start.pattern$y, window=window)
	#
	## generate the R-field
	#
	if(is.null(pre.r.field)) {
		#
		## not cyclic (toroidal) field by default
		#
		if(missing(cyclic))cyclic<-0
		#
		## lets simulate
		#
		r.field<-simulateGMRF(ncol=ncol, nrow=nrow, 
					tau=r.field.pars[2], range=r.field.pars[3], 
					nu=r.field.pars[4], expand=expand,
					xlim=window$xrange, ylim=window$yrange, cyclic=cyclic
					)
		#
		## make mean shift and exp transformation, also store parameters
		#
		r.field$v<-exp(r.field$v + r.field.pars[1])
		r.field$parameters<-r.field.pars
	}
	else{
		# 
		## check
		#
		if(!"rfhcField"%in%class(pre.r.field)) stop("pre.r.field not of class rfhcField")
		r.field<-pre.r.field
		r.field.pars<-r.field$parameters
	}
	#
	## R field ready, lets simulate points
	#
	## tweak pp for c-side compatibility
	xy<-rfhc_modify_pp(xy)
	xy$window$x<-as.numeric(r.field$xrange)
	xy$window$y<-as.numeric(r.field$yrange)
	#
	## check type of interaction
	#
	if(type=="HC")typei<-1
	else if(type=="HS") typei<-2
	else stop("type should be 'HS' (hard spheres) or 'HC' (hard core).")
	#
	## run MH simulation in c:
	#
	res<-.External("simRFHC_c", as.integer(dbg), 
			xy, r.field, 
			c(n, giveup, iter, typei),
			PACKAGE="rfhc")
	if(dbg)cat("\n")
	#
	# compile results
	#
	win2<-owin(r.field$xrange, r.field$yrange)
	z<-ppp(res[[1]], res[[2]], marks=res[[3]], window=win2)
	z$time<-format(Sys.time()-T0)
	z$r.field<-r.field
	z$parameters<-list(r.field.pars=r.field.pars, 
			r.field.model="Matern 2D approximation", 
			interaction.type=type, 
			simulation.type="MH", 
			iterations=iter)
	z<-rfhc_class_Sim(z)
	#
	##done!
	#
	if(dbg>1)print(z$time)
	#
	z
}
##
#################################################################
## modify the pp for c side compatibility
rfhc_modify_pp<-function(pp){
	n<-length(pp[["x"]])
	if(length(pp[["mass"]]) < n ) # set the masses
	{
		if(length(pp[["marks"]])< n | !is.numeric(pp[["marks"]])) pp$mass<-rep(1.0,n)
		else pp$mass<-pp$marks
	}
	if(length(pp[["types"]]) < n) # set the types
	{
		if( (is.factor(pp$marks) | is.integer(pp$marks)) & length(pp[["marks"]])==n ) pp$types<-pp$marks 
		else pp$types<-rep(1,n)
	}
	pp$mass<-as.numeric(pp$mass)
	if(length(pp$mass2)<n) pp$mass2<-rep(0.0, n)
	pp$types<-as.integer(pp$types)
	pp$marks<-NULL
	pp$window$x<-as.numeric(pp$window$x)
	pp$window$y<-as.numeric(pp$window$y)
	# include area: TODO: only rectangular?
#	pp$area<-as.numeric(area.owin(pp$window))
	pp$area<-diff(pp$window$x)*diff(pp$window$y)
	pp$x<-as.numeric(pp$x)
	pp$y<-as.numeric(pp$y)
	pp
	
	
}








