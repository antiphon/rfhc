# Fit the varying range hard-core model using Bayesian MCMC 
# 
# This one has the random jitter for Ro
#
# Author: antiphon
###############################################################################


fitRFHC_2<-function(x, iter=100, ncol=50, 
		type="HC",
		nu=2,
		mupars=c(free=1, init=NA, m=0, s2=100),
		taupars=c(free=0, init=5, F=1.05, alpha=0.25, beta=0.005),
		rangepars=c(free=0, init=0.2, F=1.05, alpha=0.25, beta=0.005),
		Rpars=list(choose=3, step=0.01, from=NULL),
		beta=9, # the Beta-likelihood nnd ~ R parameter (alfa=1)
		usePL=TRUE, # include the PL bit?
		expand, cyclic, # grid options
		seed, dbg=2, RtoNND=0.95, print.prefix="",# other options
		simulate=0) # if simulate>0, each 'simulate' step we simulate new pp    
{
	#
	## start the clock
	#
	T0<-Sys.time()
	#
	## debug printing or not
	#
	print3<-print2<-print1<-function(x){}
	if(dbg>0) print1<-function(x)cat( x, sep="")
	if(dbg>1) print2<-function(x)cat( x, sep="")
	if(dbg>2) print3<-function(x)cat( x, sep="")
	if(print.prefix!="")print0<-function()cat(print.prefix)
	else print0<-function(){}
	#
	## check type
	#
	typei<-pmatch(type, c("HC", "HS"))
	if(is.na(typei))stop("type should be 'HS'(hard spheres) or 'HC' (hard core).")
	## set functions
	#
	FUNS<-list(l=rfhc_lPL_NoPL_2, 
			   checkInteraction=rfhc_interaction_HC)
	if(usePL) FUNS$l<-rfhc_lPL_WithPL_2
	if(typei==2)FUNS$checkInteraction<-rfhc_interaction_HS
	
	#
	## aspect ratio 1
	#
	nrow<-round(ncol*diff(x$window$yrange)/diff(x$window$xrange))
	N<-nrow*ncol
	#
	## if expand missing, use template size
	#
	if(missing(expand)) expand<-nu*2+1
	Nex<-(nrow+2*expand)*(ncol+2*expand)
	#
	## seed
	#
	if(missing(seed)){runif(1); seed<-as.integer(format(Sys.time(),"%s"))*runif(1)}
	set.seed(seed)
	#		
	#
	## precompute the nearest neighbour map. This also gives us the 
	## xcol and yrow grid center points.
	#
	nndmap<-distmap(x, dimyx=c(nrow,ncol))
	xcol<-nndmap$xcol
	yrow<-nndmap$yrow
	#
	## precompute some constant stuff
	#
	CONSTS<-list(dataholders=rfhc_xy2node(x$x, x$y, xcol, yrow),
				 gridcentersxy=expand.grid(y=yrow, x=xcol),
				 x.nnd=nndist(x),
				 nndmapvalues=as.numeric(t(nndmap$v)),
				 N=N, Nex=Nex, Ndata=x$n,
				 x=x
		         ) # columnwise)
	if(missing(cyclic))cyclic<-0
	#
	## we must make sure that if we aggerage then we drop some x
	#
	# TODO	
	#
	print1("precomputations done.\n")
	#
	## check if we have too low grid
	#
	if(length(unique(CONSTS$dataholders))<x$n)
		stop("Increase resolution, more than one point hit a single pixel.")
	CONSTS$dataholders<-unique(CONSTS$dataholders)
	CONSTS$xok<-1:x$n
	if(is.null(Rpars$from))Rupdateset<-CONSTS$dataholders
	else Rupdateset<-Rpars$from
	#
	print2("global stuff set.\n")
	#
	## initial values
	#
	if(is.na(mupars[2]))mupars[2]<-log(mean(CONSTS$x.nnd*RtoNND))
	mucurrent<-mupars[2]
	taucurrent<-taupars[2]
	rangecurrent<-rangepars[2]
	#
	## initialise R
	#
	lRcurrent<-rep(0, N)
	lRcurrent[CONSTS$dataholders]<-log(CONSTS$x.nnd[CONSTS$xok]*RtoNND)
	#
	## build the state bubble for simulating log(R)
	#
	xcond<-lRcurrent[CONSTS$dataholders]-mucurrent
	rfhc_STATE<<-list(grid=list(nrow=nrow, ncol=ncol, expand=expand),
			datadists=pairdist(x), #as.matrix(dist(CONSTS$gridcentersxy[CONSTS$dataholders, ], upper=TRUE)),
			fixed=CONSTS$dataholders,
			x.cond=xcond,
			tau=taucurrent, nu=nu, range=rangecurrent,
			cyclic=cyclic, xlim=x$window$x, ylim=x$window$y,
			N=N, Nex=Nex, Ndata=length(CONSTS$dataholders),
			new=TRUE, cond.is.new=TRUE, dbg=dbg,
			IO=rfhc_insideoutside(nrow, ncol, expand)
	)
	# sample 
	rfhc_sample_conditional_GMRF() #ok 
	lRcurrent<-rfhc_STATE$x+mucurrent
	Rcurrent<-exp(lRcurrent)
	# 
	## initial Q for data
	rfhc_makeQtrue()
	#
	print2("initial values set.\n")
	#
	## log-likelihood for data
	#
	logLikData<-FUNS$l(x, Rcurrent, beta, CONSTS, FUNS)
	print2(paste("Data log start:",logLikData,"\n"))
	#
	## history 
	#
	muhist<-tauhist<-rangehist<-NULL
	pphist<-list()
	Rhist<-Rcurrent
	accRo<-accRs<-0
	#
####################  INITIAL STUFF READY LETS BEGIN! ###############################################
	for(i in 1:iter) {
		print0()
		print1(paste(i,"/", iter,": ", sep=""))
		#
		## update mu using Gibbs sampler: Normal prior -> normal posterior
		#
		if(mupars[1]){
			muA<-taucurrent*sum(rfhc_STATE$Q)+mupars[3]
			muB<-taucurrent*sum(colSums(rfhc_STATE$Q)*(lRcurrent))+mupars[3]*mupars[4]
			munew<-rnorm(1)*sqrt(1/muA)+muB/muA
			muhist[i]<-munew
		}else munew<-mucurrent
		#
		## hyper prior contribution to MH in R-udate
		#
		hyper<-0
		## update tau
		#
		if(taupars[1]){
			f<-rfhc_rlambert(taupars[3])
			taunew <- taucurrent*f
			hyper <- hyper + (taupars[4]-1)*log(f)-(taunew-taucurrent)*taupars[5]
			rfhc_STATE$cond.is.new<<-1 # because correction depends on tau
			## no need to recompute Q or L
		}else taunew<-taucurrent
		#
		## update range
		#
		if(rangepars[1]){
			f<-rlambert(rangepars[3])
			rangenew<-rangecurrent*f
			hyper <- hyper + (rangepars[4]-1)*log(f)-(rangenew-rangecurrent)*rangepars[5]
			rfhc_STATE$new<<-1   ## need to recompute Q and L
		}else rangenew<-rangecurrent
		#
		################# update R values ###################
		#
		## update R = (Ro, Rs) where Ro is the data locations and Rs others.
		#
		# First Ro:
		#
		## sample some sites to update
		## move those about
		#
		lRnewo<-lRcurrent
		if(Rpars$choose>0){
			chosen<-sample(Rupdateset, Rpars$choose)
			lRnewo[chosen]<-lRnewo[chosen]+runif(Rpars$choose, -1,1)*Rpars$step
		}
		Rnewo<-exp(lRnewo)
		#
		## check with data:
		#
		logLikDatanew<-FUNS$l(x, Rnewo, beta, CONSTS, FUNS)
		#
		## if acceptable:
		MHo<-min(exp(logLikDatanew-logLikData), 1)
		if(runif(1)<MHo) {
			logLikData<-logLikDatanew
			lRcurrent<-lRnewo
			accRo<-accRo+1
			rfhc_STATE$cond.is.new<<-1 # Ax changed
		}else {	
		}
		# Ro updated.
		### Then, given Ro 
		## sample new Rs
		rfhc_STATE$x.cond<<-lRcurrent[rfhc_STATE$fixed]-munew
		rfhc_STATE$tauold<<-taucurrent
		rfhc_STATE$tau<<-taunew
		rfhc_STATE$rangeold<<-rangecurrent
		rfhc_STATE$range<<-rangenew
		#
		rfhc_blockGMRF()
		#
		lRnew<-rfhc_STATE$x+munew
		lMHRs<-rfhc_STATE$lMH
		#
		## add hyper priors
		#
		lMHs <- lMHRs + hyper
		MHs<-min( exp(lMHs), 1 )
		#
		## now, let's see if we accept the Rs new values:
		#
		if(runif(1) < MHs) {
			lRcurrent<-lRnew
			if(taupars[1]) taucurrent<-taunew
			if(rangepars[1]) rangecurrent<-rangenew
			accRs<-accRs+1
		}
		#
		mucurrent<-munew
		Rcurrent<-exp(lRcurrent)
		#
		## store 
		#
		Rhist<-rbind(Rhist, Rcurrent)
		if(taupars[1])tauhist[i]<-taucurrent
		if(rangepars[1])rangehist[i]<-rangecurrent
		print1(paste("acc: o",accRo, " s",accRs," [MH: o", format(MHo, digits=2, width=4, zero.print=T)," s",format(MHs, digits=2, width=4, zero.print=T),"]", sep=""))
		if(taupars[1])print1(paste("[tau:",round(taucurrent,3),"]",sep=""))
		if(rangepars[1])print1(paste("[range:",round(rangecurrent,3),"]",sep=""))
		print2("\n")
		if(dbg<2)cat("             \r")
	}
	
	######################################################
	## Done. Compile results:
	#
	## results
	#
	i<-0
	mu<-mupars[2]
	if(mupars[1])mu<-c(mu,muhist)
	#	
	tau<-taupars[2]
	if(taupars[1])tau<-c(tau, tauhist)
	#
	range<-rangepars[2]
	if(rangepars[1])range<-c(range, rangehist)
	#
	opt<-data.frame(dbg=dbg, iter=iter, ncol=ncol, nrow=nrow, 
			expand=expand, seed=seed, type=type, cyclic=cyclic, usePL=usePL,
			simulate=simulate)
	took<-format(Sys.time()-T0)
	#
	results<-list(mu=mu, tau=tau, range=range, R=Rhist, time=took, 
			pphist=pphist,
			options=opt, dataholders=CONSTS$dataholders
	)
	results$acceptance<-list(Robs=accRo, Rother=accRs)
	results$parameters<-list(nu=nu, mu=mupars, tau=taupars, range=rangepars, R=Rpars, beta=beta)
	results$data<-x
	results$Rlast<-im(matrix(Rcurrent, byrow=FALSE, ncol=ncol), xcol, yrow)
	if(dbg>0)cat("Took",took,"\n")
	results<-rfhc_class_Fit(results)
	## collect options
	#
	## All done.
	results
}
## END of fit function
###########################################################
## log-pseudolikelihood function(s)
## only the independent Beta-product
#
rfhc_lPL_NoPL_2<-function(x, R, beta, CO, FU){
	if(!FU$checkInteraction(x, R, CO)) return(-Inf)
	s<-sum(log(beta)+(beta-1)*log(1-(CO$x.nnd[CO$xok]-R[CO$dataholders])/CO$x.nnd[CO$xok]))
	s
}
#
## with the Pseudo bit 
#
rfhc_lPL_WithPL_2<-function(x, R, beta, CO, FU){
	l<-rfhc_lPL_NoPL_2(x, R, beta, CO, FU)
	if(l== -Inf) return(l)
	#
	## papangelou integral (the sum part ==0 for all ok R)
	# 
	Rok<-R<CO$nndmapvalues
	ok<-apply(cbind(CO$gridcentersxy[Rok,], R[Rok]), 1, function(xy) sum( ((xy[1]-x$x)+(xy[2]-x$y))< xy[3]^2)<1 )
	l-area.owin(CO$x$window)*sum(ok)/CO$N
}
#############################################################
## intearction, HS or HC
rfhc_interaction_HC<-function(x, R, CON){
	a<-CON$x.nnd[CON$xok] < R[CON$dataholders]
	!(sum(a)>0)
}
rfhc_interaction_HS<-function(x, R, CON){
	print("Not implemented")
}

############################################################  
## Lambert's W distribution:
rfhc_rlambert<-function(f){
	A<-f+2*log(f)-1/f
	B<-(log(f)-1/f)/A
	x<-exp(A*(runif(1)-B))
	eps<-w<-1  
	while(eps>0.001)eps<-abs(w-(w<-(x*exp(-w)+w^2)/(w+1)))
	w
}
##
