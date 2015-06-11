# Trickery with grid
# 
# Author: tarajala
###############################################################################

####################################################
## take nrow, ncol, expand which represent an 
## (nrow+2*expand)x(ncol +2*expand) grid.
##
## Turn the (i,j) grid to nodes.
##
## Return a list with two components:
## IN = the inside node indices
## OUT = the ... yeah
##
rfhc_insideoutside<-function(nrow, ncol, expand) {
	in_ij<-expand.grid(1:nrow+expand, 1:ncol+expand)
	nyex<-nrow+2*expand
	nxex<-ncol+2*expand
	IN<-apply(in_ij, 1, function(ij)rfhc_grid2node(ij[1], ij[2], nyex, nxex))
	OUT<-setdiff(1:(nyex*nxex), IN)
	list(IN=IN, OUT=OUT)
}
## 
###################################################
###################################################
## The grid to node to grid stuff
###################################################
rfhc_node2grid<-function(k, nrow, ncol){
	j<-trunc((k-1)/nrow)
	c(k-j*nrow, j+1)
}
###################################################
rfhc_grid2node<-function(i, j, nrow, ncol){
	i+(j-1)*nrow	
}
###################################################
## the xy to grid/node stuff.
###################################################
rfhc_xy2grid<-function(x, y, xcol, yrow){
	dx<-diff(xcol[1:2])/2
	dy<-diff(yrow[1:2])/2
	i<-apply(cbind(y+dy), 1, function(yi) sum(yi>yrow) )
	j<-apply(cbind(x+dx), 1, function(xi) sum(xi>xcol) )
	cbind(i, j)
}
###################################################
rfhc_xy2node<-function(x, y, xcol, yrow){
	ij<-rfhc_xy2grid(x, y, xcol, yrow)
	rfhc_grid2node(ij[,1], ij[,2], length(yrow), length(xcol))
}
###################################################

