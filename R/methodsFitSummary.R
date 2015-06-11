plot.rfhcFitSummary<-function(x, ..., maxcols=4) {
	pred<-!is.null(x$Rpred$mean.im)
	n <- x$parameters$mu[1]+x$parameters$tau[1]+x$parameters$range[1] + 2*pred
	ncol<-min(maxcols, n)
	nrow<-ceiling(n/ncol)
	par(mfrow=c(nrow,ncol))
	pars<-x$parameters
	if(!is.null(x$mu)){
		plot(x$mu$density, main="Mean of the log-field")
		curve(dnorm(x, pars$mu[3], pars$mu[4]), add=TRUE, lty=2, col=2)
	}
	if(!is.null(x$tau)){
		plot(x$tau$density, main="Precision of the log-field")
		curve(dgamma(x, pars$tau[4], pars$tau[5]), add=TRUE, lty=2, col=2)
	}
	if(!is.null(x$range)){
		plot(x$range$density, main="Range of the log-field")
		curve(dgamma(x, pars$range[4], pars$range[5]), add=TRUE, lty=2, col=2)
	}
	if(pred){
		plot(x$Rpred$mean.im, main=paste("Posterior predictive mean"), ...)
		plot(x$Rpred$sd.im, main="Posterior predictive sd", ...)
	}
}

print.rfhcFitSummary<-function(x, ...) {
	cat("Summary for rfhc-fit:\n")
	om<-setdiff(names(x$options),c("predict.x","predict.R","type","usePL","dbg"))
	cat("* Options | ", paste(names(x$options[om]),x$options[om],sep=":"),"\n")
	cat("* Predictions every nth iteration | x:", x$options$predict.x, " R:", x$options$predict.R,"\n",sep="")
	cat("* Chain sampling | ", paste(names(x$sampling),x$sampling,sep=":"),"")
	cat("=> (hyper)n=", floor((x$options$iter-x$sampling[1]+1)/x$sampling[2])[[1]],sep="")
	cat(" (R)n=",floor((floor(x$options$iter/x$options$predict.R)-x$sampling[3]+1)/x$sampling[4])[[1]],sep="","\n")
	cat("* R update | ")
	cat("choose:",x$param$R$choose,sep="")
	cat(" from:")
	if(is.null(x$param$R$from))cat("all\n")
	else cat("given subset\n")
	cat("* Beta-parameter | ", x$param$beta,"\n")
	cat("* Hyper parameters:\n")
	cat(" * Fixed | ")
	cat("nu:",x$param$nu)
	pars<-c("mu","tau","range")
	for(p in pars)if(x$param[[p]][1]==0)cat(" ",p, ": ", x$param[[p]][2],sep="")
	cat("\n * Estimated:\n")
	est<-NULL
	for(e in x$estimated)est<-rbind(est, x[[e]]$summary)
	rownames(est)<-x$estimated
	print(est)
	
	
}