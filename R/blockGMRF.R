# Block update for GMRF + hypers
# 
# Author: antiphon
###############################################################################

rfhc_blockGMRF<-function() {
	#
	## range in real scale
	## mean 0
	#
	# we use algoritm 2.6. in rue&held 05
	#
	# this function simulates new x and computes lMH
	#
	## are Q and L computed anew?
	is.new<-rfhc_STATE$new
	#
	## compute Q for the data if range has changed
	#
	rfhc_STATE$Qdata.old<<-rfhc_STATE$Qdata
	if( is.new ) {
		rfhc_makeQtrue()
	}
	#
	## sample new x if we want to predict
	#
	if(rfhc_STATE$predict | is.new) rfhc_sample_conditional_GMRF()
    if((rfhc_STATE$predict | is.new) & rfhc_STATE$use.RandomFields) rfhc_sample_conditional_GRF()
	#
	## compute the lMH, depends only on data
	#
	# if range changed
	if(is.new) {
		a<-rfhc_STATE$Ndata*log(rfhc_STATE$tau/rfhc_STATE$tauold)/2 + 0.5*(rfhc_LOGDET(rfhc_STATE$Qdata)-rfhc_LOGDET(rfhc_STATE$Qdata.old))
		b<- -0.5*rfhc_STATE$tau*t(rfhc_STATE$x.cond)%*%rfhc_STATE$Qdata%*%rfhc_STATE$x.cond + 0.5*rfhc_STATE$tauold*t(rfhc_STATE$x.cond)%*%rfhc_STATE$Qdata.old%*%rfhc_STATE$x.cond
	}
	else {
		a<- rfhc_STATE$Ndata*log(rfhc_STATE$tau/rfhc_STATE$tauold)/2
		b<- -0.5*(rfhc_STATE$tau-rfhc_STATE$tauold)*t(rfhc_STATE$x.cond)%*%rfhc_STATE$Qdata%*%rfhc_STATE$x.cond
	}
	rfhc_STATE$lMH<<- as.numeric(a + b)		
}
#####################################################################3
# didnt work
#	else { 
#		a<-(rfhc_STATE$N-rfhc_STATE$Ndata)*log(rfhc_STATE$tau/rfhc_STATE$tauold)		
#		o<-rfhc_STATE$others
#		f<-rfhc_STATE$fixed
#		S <-      t(rfhc_STATE$x.old) %*% rfhc_STATE$Q %*% rfhc_STATE$x.old
#		Stilde <- t(rfhc_STATE$x)     %*% rfhc_STATE$Q %*% rfhc_STATE$x
#		
#		Supp <-t(c(rfhc_STATE$x.old[o],rfhc_STATE$x[f]))  %*% rfhc_STATE$Q[c(o,f),c(o,f)] %*% c(rfhc_STATE$x.old[o],rfhc_STATE$x[f])		
#		Sdown <-t(c(rfhc_STATE$x[o],rfhc_STATE$x.old[f])) %*% rfhc_STATE$Q[c(o,f),c(o,f)] %*% c(rfhc_STATE$x[o],rfhc_STATE$x.old[f])
#		a1<- -0.5*rfhc_STATE$tau*(Stilde + Supp)
#		b1<-  0.5*rfhc_STATE$tauold*(S + Sdown)
#		
#		a2<- -rfhc_STATE$tauold * t(rfhc_STATE$x.old[f])%*% rfhc_STATE$Qdata%*%rfhc_STATE$x.old[f]
#		b2<- rfhc_STATE$tau*t(rfhc_STATE$x[f])%*%rfhc_STATE$Qdata%*%rfhc_STATE$x[f]
#		rfhc_STATE$lMH<<- as.numeric(a + a1 + a2 + b1 + b2)
#	}
#}



############################################################
## MH is computed only the Rs update
#rfhc_blockGMRF_just_Rs<-function() {
#	#
#	## range in real scale
#	## mean 0
#	#
#	# we use algoritm 2.6. in rue&held 05
#	#
#	# this function simulates new x and computes lMH
#	#
#	## are Q and L computed anew?
#	is.new<-rfhc_STATE$new
#	#
#	## compute Q for the data if range has changed
#	#
#	if( is.new ) {
#		rfhc_STATE$Q.old<<-rfhc_STATE$Q.data
#		rfhc_makeQtrue()
#		print("[Qdata]")
#	}
#	#
#	## sample new x
#	#
#	rfhc_sample_conditional_GMRF()
#	#
#	## compute the lMH
#	#
#	# if range changed
#	if(is.new) {
#		print("RANGE CHANGE NOT IMPLEMENTED")
#	}
#	else {
#	a<-log( rfhc_STATE$tau/rfhc_STATE$tauold )*rfhc_STATE$Ndata*0.5
#	b<- - 0.5 * (rfhc_STATE$tau-rfhc_STATE$tauold) * 
#			as.numeric( t(rfhc_STATE$x.cond)%*%rfhc_STATE$Qdata%*%rfhc_STATE$x.cond )
#	rfhc_STATE$lMH<<- a + b
#    }
#}



