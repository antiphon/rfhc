# Methods for class rfhcFit
# 
# Author: antiphon
###############################################################################

plot.rfhcFit<-function(x, burnin=1, thin=1, Rburnin=1, Rthin=1, data=TRUE, sigma=FALSE, ...) {
	#
	## check if we're dealing with a simulation data
	## plot extra diagnostics.
	# 
	#
	sim<-"rfhcSim"%in%class(x$data)
	if(sim) mmapOK<-require(mmap, quietly=TRUE)
	else mmapOK<-FALSE
	# sampling
	these<-seq(burnin, x$options$iter, by=thin)
	Rthese<-seq(Rburnin, nrow(x$Rpred), by=Rthin)
	if(length(Rthese)<3)pred<-FALSE
	else pred<-TRUE
	#
	## layout
	#
	n<-x$parameters$mu[1]+x$parameters$tau[1]+x$parameters$range[1]+1 + 1*data + 1*sim*data
	nrow<-ceiling(n/3)
	ncol<- ifelse(nrow==1,n,3)
	#
	if(n>1)par(mfrow=c(nrow, ncol))
	#
	## if variance instead of precision
	#
	if(sigma) x$tau<-1/x$tau
	#
	# plot data 
	#
	pp2<-pp<-x$data
	Ro<-x$Robs[these,]
	v<-apply(Ro, 2, median)
	
	if(data){
		if(sim)plot(pp, main="data", size=1, ...)
		else plot(pp, main="data")
		#
		## plot squared difference plot if we know the field
		# 
		if(sim & mmapOK){
			SS<-((pp$marks-v)^2)
			pp0<-pp
			pp0$marks<-SS
			SSmap<-meanmap(pp0, adaptive=TRUE, kernel=2)
			plot(SSmap, main="squared difference from median R to true R")
			symbols(circles=SS, pp$x, pp$y, inches=FALSE, add=TRUE)
		}
	}
	## plot median R mark pattern 
	## plot on top of predictive mean field if applicable 
	#
	pp2$marks<-v
	tit<-"Median marks"
	if(pred) {
		Rpred<-x$Rlast
		Rpred$v<-matrix(c(apply(x$Rpred[Rthese,], 2, median)), ncol=ncol(Rpred$v))
		plot(Rpred, ...)
		tit<-"Posterior predictive field + median marks"
	}
	plot.ppp(pp2, use.marks=TRUE, markscale=1, add=pred, inches=FALSE, main="")
	points(pp2$x,pp2$y, cex=0.2, pch=19)
	title(tit)
	## plot chains for the others
	for(b in c("mu","tau","range")) {
		if(x$parameters[[b]][1]>0) {
			ts.plot(x[[b]][these], main=b, ylab="")
			if(sim) abline(h=pp$parameters$r.field.pars[which(b==c("mu","tau","range"))], lty=2, col=gray(.4))
		}
	}
	
}

#########################################################################
plotRtraj.rfhcFit<-function(x, these, burnin=1, thin=1, ...) {
	if(missing(these)) these<-1:ncol(x$Robs)
	samp<-seq(burnin, nrow(x$Robs), by=thin)
	plot(NA, xlab="iter", ylab="R_i", xlim=c(0,length(samp)), ylim=c(0,max(x$Robs[,these])))
	for(j in these)lines(x$Robs[samp,j])
}

#########################################################################
plotRpost.rfhcFit<-function(x, burnin=1, thin=1, these=1:6, maxcol=6, breaks=50, error=TRUE, ...) {
	## plot the posterior of R-dataholders
	samp<-seq(burnin, x$options$iter, by=thin)
	den<-function(z) hist(z, breaks=breaks, col="black",  freq=FALSE, main="")
	
	n<-length(these)
	
	par(mfrow=c(ceiling(n/maxcol), min(n,maxcol)), mar=c(2,1,1,1), cex=0.4)
	
	beta<-x$parameters$beta
	
	nnd<-nndist(x$data)
	#
	for(i in these){
		j<-i
		nd<-nnd[j]
		if(error){
			v<-(nd-x$Robs[samp,i])/nd
			den(v)
			curve(dbeta(x, 1, beta), add=TRUE, col=2)
		}
		else{
			v<-x$Robs[samp,i]
			den(v)
			rug(nd, lwd=2, col=3)
		}
		
	}
}

#########################################################################
#########################################################################

summary.rfhcFit<-function(x, burnin=1, thin=1, Rburnin=1, Rthin=1, paradj=1, Radj=1, ...) {
	these<-seq(burnin, x$options$iter, by=thin)
	Rthese<-seq(Rburnin, nrow(x$Rpred), by=Rthin)
	pred<-ifelse(length(Rthese)>2, TRUE, FALSE)
	result<-list()
	result$parameters<-x$parameters
	result$options<-x$options
	result$sampling<-data.frame(burnin=burnin, thin=thin, Rburnin=Rburnin, Rthin=Rthin)
	## helper functions
	qfun<-function(z)quantile(z, probs=c(0.05,0.95))
	sumup<-function(z) data.frame(mean=mean(z), median=median(z), 
				sd=sd(z), q=rbind(qfun(z)), min=min(z),max=max(z))
	den<-function(z) density(z, adjust=paradj, n=125)
	denR<-function(z){zz<-density(z, adjust=Radj, n=125); cbind(x=zz$x, y=zz$y)}
	##
	result$estimated<-NULL
	## parameters
	for(p in c("mu","tau","range")){
		if(x$parameters[[p]][1]) {
			v<-list()
			v$summary<-sumup(x[[p]][these])
			v$density<-den(x[[p]][these])
			result[[p]]<-v
			result$estimated<-c(result$estimated, p)
		}
	}
	## R obs:
	Robs<-Rpred<-list()
	## density
	v<-lapply(split(t(x$Robs[these, ]), factor(1:length(x$dataholders))), denR)
	names(v)<-as.character(x$dataholders)
	#
	## obs summary
	#
	Robs$density<-v
	v<-apply(x$Robs[these, ], 2, sumup)
	Robs$summary<-t(sapply(v, c))
	rownames(Robs$summary) <- as.character(x$dataholders)
	#
	## predictive mean and variance maps of Rpred
	#
	if(pred){
		v<-apply(x$Rpred[Rthese, ], 2, mean)
		Rpred$mean.im<-x$Rlast
		Rpred$mean.im$v<-matrix(v, ncol=ncol(Rpred$mean.im))
		v<-apply(x$Rpred[Rthese, ], 2, sd)
		Rpred$sd.im<-x$Rlast
		Rpred$sd.im$v<-matrix(v, ncol=ncol(Rpred$sd.im))
		#
		## variance estimate
		#
		v<-c(apply(x$Rpred[Rthese, ], 2, function(x)var(log(x))))
		Rpred$var<-list(density=denR(v), summary=sumup(v))
		
	}
	result$Robs<-Robs
	result$Rpred<-Rpred
	## mean square difference if data is simulation:
	if("rfhcSim"%in%class(x$data)){
		result$MSE_mean <- mean((unlist(Robs$summary[,"mean"])-x$data$marks)^2)
		result$MSE_median <- mean((unlist(Robs$summary[,"median"])-x$data$marks)^2)
	}
	rfhc_class_FitSummary(result)	
}
#########################################################################
rfhc_fit_to_im<-function(x, burnin=1, thin=1, fun=mean, ...){
	these<-seq(burnin, x$options$iter, by=thin)
	R<-matrix(apply(x$R[these,], 2, fun, ...), byrow=!TRUE, ncol=x$options$ncol)
	w<-x$data$window
	dx2<-diff(w$xrange)/(2*x$options$ncol)
	dy2<-diff(w$yrange)/(2*x$options$nrow)
	xcol<-seq(w$xrange[1]+dx2, w$xrange[2]-dx2, length=x$options$ncol)
	yrow<-seq(w$yrange[1]+dy2, w$yrange[2]-dy2, length=x$options$nrow)
	rfhc_class_Field(im(R, xcol=xcol, yrow=yrow, xrange=x$data$window$xrange, yrange=x$data$window$yrange))
	
}
####################################################
## compute the posterior mean/median at 
## data points with post. means of parameters 
## and predict elsewhere
#predict.rfhcFit<-function(x, burnin=1, thin=1, fun=mean, ...) {
#	these<-seq(burnin, x$options$iter, by=thin)
#	p<-x$parameters
#	# 
#	## collect field parameters
#	if(p$mu[1]) mu<-fun(x$mu[these])
#	else mu<-p$mu[2]
#	if(p$tau[1]) tau<-fun(x$tau[these])
#	else tau<-p$tau[2]
#	if(p$range[1]) range<-fun(x$range[these])
#	else range<-p$range[2]
#	nu<-p$nu
#	#
#	## post mean/median values at x
#	#
#	v<-log(apply(x$Robs[these, ], 2, fun)-mu)
#	#
#	cat("Simulation with tau=",tau," range=",range, " nu=",nu, "\n",sep="")
#	#
#	z<-simulate_condGMRF(tau=tau, range=range, nu=nu, 
#			xlim=x$data$window$x, ylim=x$data$window$y,
#			data=list(x=x$data$y, y=x$data$y, v=v), ...)
#	z$v<-exp(z$v+mu)
#	z
#}
##########################################################
######## Model fit
##
## error maps
## x = fit
##
error.mapfun<-function(x, burnin=1, thin=1, mapfun=distmap, ...) {
	dmaps<-list()
	# sample
	these<-seq(burnin, length(x$pphist), by=thin)
	if(length(these)<5) warning("Less than 5 datapoints.")
	# compute
	k<-0
	for(i in these)dmaps[[k<-k+1]]<-mapfun(x$pphist[[i]], ...)
	dmapd<-mapfun(x$data,  ...)
	ncol<-dmapd$dim[2]; nrow<-dmapd$dim[1]
	# scaled mean square errors per pixel
	dv<-c(dmapd$v)
	v<-sapply(dmaps, function(x)(c(x$v)-dv)^2)
	error.MMSE<-error.MSE<-error.MSE_std<-error.sd<-dmapd
	error.MSE$v<-matrix(c(apply(v, 1, mean)), ncol=ncol)
	error.MSE_std$v<-matrix(c(apply(v, 1, function(x)mean(x)/sd(x))), ncol=ncol)
	error.sd$v<-matrix(c(apply(v, 1, sd)), ncol=ncol)
	# compute the Integrated MSE i.e. sum over the pixelwise MSE:s
	# and MMSE which is mean pixel MSE.
	vm<-apply(sapply(dmaps, function(x)x$v), 1, mean)
	vm_sd<-apply(sapply(dmaps, function(x)x$v), 1, sd)/sqrt(length(these))
	MMSE<-mean((dv - vm)^2)
	MMSE_std<-mean((dv - vm)^2/vm_sd)
	IMSE<-sum(error.MSE$v)*area.owin(x$data$window)/(ncol*nrow)
	IMSE_std<-sum(error.MSE_std$v)*area.owin(x$data$window)/(ncol*nrow)
	cat("IMSE (std)=", IMSE,"(", IMSE_std,")\n")
	cat("MMSE (std)=", MMSE,"(", MMSE_std,")\n")
	# done.	
	list(predicted=dmaps, data=dmapd, error.MSE=error.MSE, error.MSE_std=error.MSE_std, error.sdMSE=error.sd, IMSE=IMSE, IMSE_std=IMSE_std, MMSE=MMSE, MMSE_std=MMSE_std, n=length(these))
}

## error summaries
## x = fit
## use fv-class so that we pick the mean values
##
error.sumfun<-function(x, burnin=1, thin=1, sumfun=pcf, ...) {
	dsums<-list()
	# sample
	these<-seq(burnin, length(x$pphist), by=thin)
	if(length(these)<5) warning("Less than 5 datapoints.")
	# compute
	k<-0
	for(i in these)dsums[[k<-k+1]]<-sumfun(x$pphist[[i]], ...)
	datas<-sumfun(x$data, ...)
	# done.	
	list(predicted=dsums, data=datas, n=length(these))
}

#########################################
print.rfhcFit<-function(x, ...){
	cat("Result object of fitRFHCmodel.\n")
	cat("* Data: \n\n")
	print(x$data)
	om<-setdiff(names(x$options),c("predict.x","predict.R","type","usePL","dbg"))
	cat("\n* Options | ", paste(names(x$options[om]),x$options[om],sep=":"),"\n")
	cat("* Predictions every nth iteration | x:", x$options$predict.x, " R:", x$options$predict.R,"\n",sep="")
	cat("* Length of prediction chains | x:", length(x$pphist), " R:", nrow(x$Rpred),"\n")
	cat("* R update | ")
	cat("choose:",x$param$R$choose,"(",round(x$param$R$choose/x$data$n*100),"%)", sep="")
	cat(" from:")
	if(is.null(x$param$R$from))cat("all\n")
	else cat("given subset\n")
	cat("* Beta-parameter | ", x$param$beta,"\n")
	cat("* Hyper parameters:\n")
	cat(" * Fixed | ")
	cat("nu:",x$param$nu)
	pars<-c("mu","tau","range")
	for(p in pars)if(x$param[[p]][1]==0)cat(" ",p, ": ", x$param[[p]][2],sep="")
	cat("\n * Free |")
	for(p in pars)if(x$param[[p]][1]>0)cat(" ",p,sep="")
	cat("\n* Fitting took | ",x$time,"\n")
}