library(pdist);

weighted.geno=function(geno,weight) {
	WX=matrix(,nrow=nrow(geno),ncol=ncol(geno));
	for (j in 1:ncol(WX)) {
		WX[,j]=weight[j]*geno[,j]
	}
	return(WX);
}

# h is a positive bandwidth parameter;
# h controls the rate of decay of the genetic correlation between two individuals
within.env.WK=function(weighted.geno,h) {
	ED=as.matrix( (dist(weighted.geno,method='euclidean'))^2 ) 
	ED=ED/max(ED); # re-scale similarity into range (0,1) 
	kernel=exp(-h*ED);
	return(kernel);
}

cross.env.WK=function(weighted.geno1,weighted.geno2,h) {
	ED=(as.matrix(pdist(weighted.geno1,weighted.geno2)))^2; 
	ED=ED/max(ED); # re-scale similarity into range (0,1)
	kernel=exp(-h*ED);
	return(kernel);
}
