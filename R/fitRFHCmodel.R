# Fit the varying range hard-core model using Bayesian MCMC 
#
# This one has the beta-distribution draw for Ro update
# Has also a different likelihood
#
# 
# Author: antiphon
###############################################################################


fitRFHC<-function(x, iter=100, ncol=50, 
		nu=2,
		range=0.3,
		mupars=c(free=1, init=NA, m=0, s2=100),
		taupars=c(free=1, init=5, F=1.05, alpha=0.25, beta=0.005),
		Rpars=list(choose=0.15, from=NULL), ## 15% of observed R's updated per iteration
		beta=9, # the Beta-likelihood nnd ~ R parameter (alfa=1)
		expand, cyclic, # grid options
		seed, dbg=1, RtoNND=0.999, print.prefix="", type="HC",# other options
		predict.R=10, # predict R everywhere every nth-step.
		predict.x=0  # if predict.x>0, each 'simulate' step we simulate new pp
	    )#use.RandomFields=FALSE)
{
	use.RandomFields <- FALSE # set true if you want to try out the randomFields package for conditional simulation.
	#
	## PL is not working since we dont have the R outside data. So set the following to FALSE 
	#
	usePL<-FALSE # include the PL bit?
	#
	## range should be fixed...
	#
	rangepars<-range
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
	################################################
	# check parameters
	#
	## check type
	#
	typei<-pmatch(type, c("HC", "HS"))
	if(is.na(typei))stop("type should be 'HS'(hard spheres) or 'HC' (hard core).")
	## set functions
	#
	FUNS<-list(l=rfhc_lPL_NoPL, 
			   checkInteraction=rfhc_interaction_HC)
	if(usePL) FUNS$l<-rfhc_lPL_WithPL
	if(typei==2)FUNS$checkInteraction<-rfhc_interaction_HS
	#
	## field parameters
	#
	## handle single fixed values
	if(length(mupars)==1) mupars<-c(0,mupars)
	if(length(taupars)==1) taupars<-c(0, taupars)
	if(length(rangepars)==1) rangepars<-c(0, rangepars)
	## handle free but only intialised (no priors)
	if(length(mupars)==2 & mupars[1]) mupars<-c(1, mupars[2], m=0, s2=100)
	if(length(taupars)==2 & taupars[1]) taupars<-c(1, taupars[2], F=1.1, alpha=0.25, beta=.005)
	if(length(rangepars)==2 & rangepars[1]) rangepars<-c(1, rangepars[2], F=1.05, alpha=0.25, beta=.005)
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
	
	## TEMP
#	step<-1/ncol
#	xcolex<-yrowex<-seq(0-expand*step,1+expand*step,length=ncol+2*expand)
#	extendedgriddists<<- as.matrix(dist( expand.grid(y=yrowex, x=xcolex)  ,upper=T))
	#
	## seed
	#
	if(missing(seed)){runif(1); seed<-as.integer(format(Sys.time(),"%s"))*runif(1)}
	set.seed(seed)
	#
	## check predictive paramaters
	#
	if(predict.x>0 & predict.x != predict.R){
		warning("Point pattern predictive simulation requested but R-field prediction step does not match. Additional R-field predictions will be performed.")
	}
	############################################################
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
	#
	if(is.null(Rpars$from))Rupdateset<-CONSTS$dataholders
	else Rupdateset<-Rpars$from
	#
	if(is.na(Rpars$choose))Rpars$choose<-CONSTS$Ndata
	# percentual choice
	if(Rpars$choose<1 & Rpars$choose>0) Rpars$choose<-round(CONSTS$Ndata*Rpars$choose)
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
			datadists=pairdist(x), 
			fixed=CONSTS$dataholders,
			others=setdiff(1:CONSTS$N, CONSTS$dataholders),
			x.cond=xcond,
			tau=taucurrent, nu=nu, range=rangecurrent,
			cyclic=cyclic, xlim=x$window$x, ylim=x$window$y,
			N=N, Nex=Nex, Ndata=length(CONSTS$dataholders),
			new=TRUE, cond.is.new=TRUE, dbg=dbg, predict=TRUE,
			IO=rfhc_insideoutside(nrow, ncol, expand),
            use.RandomFields=use.RandomFields
	     )
    ### add data if RandomFields used
    if(use.RandomFields) rfhc_STATE$data<<-x
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
	## use separate histories for data and prediction
	Robshist<-Rcurrent[rfhc_STATE$fixed]
	Rpredhist<-Rcurrent
	accRo<-accRs<-0
	#
	#
	## starting values etc print.
	#
	hyp<-paste(ifelse(mupars[1],"mu",""), ifelse(taupars[1],"tau",""), ifelse(rangepars[1],"range","") )
	print1(paste("* Start rfhc-model fitting,",date(),":\n* R(x) values updated per iteration:",Rpars$choose," (",round(100*Rpars$choose/CONSTS$Ndata, 1),"%)\n* Estimated hyper parameters:",hyp,
					"\n* Predictions at every (nth) step: x (",predict.x,") R (",predict.R,")\n***\n"))
	
####################  INITIAL STUFF READY LETS BEGIN! ###############################################
	for(i in 1:iter) {
		print0()
		print1(paste(i,"/", iter,": ", sep=""))
		## check prediction request
		#
		# set R field prediction OFF by default
		rfhc_STATE$predict<<-FALSE
		# if predictions requested
		if((predict.R>0) & (i%%predict.R==0) ) rfhc_STATE$predict<<-TRUE
		# if predictive point pattern simulation requested we must update the R-field everywhere.
		if((predict.x>0) & (i%%predict.x==0) ) rfhc_STATE$predict<<-TRUE
		#
		## update mu using Gibbs sampler: Normal prior -> normal posterior
		#
		if(mupars[1]){
			muA<-taucurrent*sum(rfhc_STATE$Qdata)+mupars[3]
			muB<-taucurrent*sum(colSums(rfhc_STATE$Qdata)*(lRcurrent[rfhc_STATE$fixed]))+mupars[3]*mupars[4]
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
			f<-rfhc_rlambert(rangepars[3])
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
		lRnew<-lRcurrent
		if(Rpars$choose>0){
			## ## TODO: clean up 
			## ## NEW: sample indiv Roi from beta(b, di)
			## ## NOTE: also the likelihood is changed!
			#
			chosen<-sample(Rupdateset, Rpars$choose)
			## this awful trickery is to keep some nodes fixed...
			notchosen  <- setdiff(rfhc_STATE$fixed, chosen)
			temp       <- lRnew[notchosen]
			u          <- rbeta(rfhc_STATE$Ndata, 1, beta)
			newones    <- log(CONSTS$x.nnd*(1-u))
			lRnew[rfhc_STATE$fixed]   <- newones
			lRnew[notchosen]          <- temp
			rfhc_STATE$cond.is.new <<- 1
		} ## new Ro's sampled
		# MH for the Ro's 
		#
		## HC + PL
		#
		Rnew<-exp(lRnew)
		logLikDatanew<-FUNS$l(x, Rnew, beta, CONSTS, FUNS)
		# check
		if(logLikDatanew==-Inf)print("Something wrong: Accepting HC violating R-field.")
		#
		## Ro's MH-bit from log-Gaussian field 
		#
		rfhc_STATE$x.cond.old <<- lRcurrent[rfhc_STATE$fixed]-munew
		rfhc_STATE$x.cond     <<- lRnew[rfhc_STATE$fixed]-munew
		rfhc_STATE$tau        <<- taucurrent
		rfhc_lMHData()
		#
		#  add HC+PL  and log-Gaussian
		#
		lMHo<-rfhc_STATE$lMHdata + logLikDatanew-logLikData + sum(lRcurrent[rfhc_STATE$fixed] - lRnew[rfhc_STATE$fixed])
		MHo <- min(exp(lMHo), 1)
		if(runif(1) < MHo) {
			accRo<-1+accRo
			lRcurrent[rfhc_STATE$fixed]<-lRnew[rfhc_STATE$fixed]
		}
		else{
			rfhc_STATE$x.cond<<-rfhc_STATE$x.cond.old
		}
		#
		### Then, given new Ro 
		## sample new Rs (if requested) and compute the log-MH for parameters theta.
		#
		## set the problem
		rfhc_STATE$tauold<<-taucurrent
		rfhc_STATE$tau<<-taunew
		rfhc_STATE$rangeold<<-rangecurrent
		rfhc_STATE$range<<-rangenew
		#
		rfhc_blockGMRF()
		#
		lRnew <- rfhc_STATE$x + munew
		# in case no update of x
		lRnew[rfhc_STATE$fixed] <- rfhc_STATE$x.cond + munew
		#
		#
		## MH: hyper priors and data log and Rs log
		#
		lMHRs <- rfhc_STATE$lMH
		lMHs   <- lMHRs + hyper
		MHs    <- min( exp(lMHs), 1 )
		#
		## now, let's see if we accept the Rs new values:
		#
		rfhc_STATE$accepted<<-0	
		if(runif(1) < MHs) {
			lRcurrent<-lRnew
			if(taupars[1]) taucurrent<-taunew
			if(rangepars[1]) rangecurrent<-rangenew
			accRs<-accRs+1
			rfhc_STATE$accepted<<-1
		}
		#
		mucurrent<-munew
		Rcurrent<-exp(lRcurrent)
		#
		#
		## store 
		#
		Robshist<-rbind(Robshist, Rcurrent[rfhc_STATE$fixed])
		if(taupars[1])tauhist[i]<-taucurrent
		if(rangepars[1])rangehist[i]<-rangecurrent
		print1(paste("acc: o",accRo, " s",accRs," [MHo: ", format(MHo, digits=2, width=4, zero.print=T), " MHs:", format(MHs, digits=2, width=4, zero.print=T), "]", sep=""))
		if(taupars[1])print1(paste("[tau:",round(taucurrent,3),"]",sep=""))
		if(rangepars[1])print1(paste("[range:",round(rangecurrent,3),"]",sep=""))
		#
		## simulate a point pattern from predictive distribution if needed
		#
		if( (predict.x>0) & (i%%predict.x==0) ) {
				print1("predict.x[")
				rim<-rfhc_class_Field(im(matrix(Rcurrent, ncol=ncol, byrow=FALSE), xcol, yrow))
				pphist<-append(pphist, list(simulateRFHC(n=x$n, pre.r.field=rim, window=x$window, type=type)))
				print1("]")
				gc()
				Rpredhist<-rbind(Rpredhist, Rcurrent)
		}
		else if((predict.R>0) & (i%%predict.R==0) ){Rpredhist<-rbind(Rpredhist, Rcurrent)}
		
		#		
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
			predict.x=predict.x, predict.R=predict.R)
	took<-format(Sys.time()-T0)
	#
	results<-list(mu=mu, tau=tau, range=range, Robs=Robshist, Rpred=rbind(Rpredhist), time=took, 
			pphist=pphist,
			options=opt, dataholders=CONSTS$dataholders
	)
	results$acceptance<-list(Robs=accRo, other=accRs)
	results$parameters<-list(nu=nu, mu=mupars, tau=taupars, range=rangepars, R=Rpars, beta=beta)
	results$data<-x
	results$Rlast<-im(matrix(Rcurrent, byrow=FALSE, ncol=ncol), xcol, yrow)
	if(dbg>0)cat("\nTook",took,"\n")
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
rfhc_lPL_NoPL<-function(x, R, beta, CO, FU){
	if(!FU$checkInteraction(x, R, CO)) return(-Inf)
	s<-0
	s
}
#
## with the Pseudo bit 
#
rfhc_lPL_WithPL<-function(x, R, beta, CO, FU){
	l<-rfhc_lPL_NoPL(x, R, beta, CO, FU)
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
#############################################################
## compute the lMH ratio for the Ro-update
rfhc_lMHData<-function() {
	a <- t(rfhc_STATE$x.cond)%*%rfhc_STATE$Qdata%*%rfhc_STATE$x.cond
	b <- t(rfhc_STATE$x.cond.old)%*%rfhc_STATE$Qdata%*%rfhc_STATE$x.cond.old
	rfhc_STATE$lMHdata<<--0.5*rfhc_STATE$tau * ( a-b )
}
