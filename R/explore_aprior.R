#########################################################
# Functions that explore data for a priori
# 
#

##################################
## variogram
rfhc.vgram<-function(x, indep=TRUE, N=20, ..., p=TRUE){
	OK<-require(fields, quietly=TRUE)
	if(!OK) stop("package 'fields' required for variogram.")
	
	
	v<-vgram(cbind(x$x, x$y), log(nndist(x)), N=N, ...)
	
	if(p)plot(v$centers, v$stats["mean",], col=1,type="b", pch=19, cex=0.3, main="2*variogram")
	
	if(p& indep & "rfhcSim" %in%class(x)) {
		xy<-expand.grid(x$r.field$xcol, x$r.field$yrow)
		
		for(i in 1:(1*indep)){
			i<-sample(1:nrow(xy), min(nrow(xy), 2*x$n))
			v0<-vgram(xy[i,], log(x$r.field$v[i]), ...)
			lines(v0$centers, v0$stats["mean",], col=2,type="b", pch=2, cex=0.3)
		}
	}
	v
}