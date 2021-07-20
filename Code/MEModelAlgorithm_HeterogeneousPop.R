library(matrixcalc); 
library(coda);
library(geoR) 
library(mvnfast) 

library(foreach)
library(doParallel);
library(parallel);
library(iterators);

## Y[j] is phenotype data of the jth environment
## K[j] is the genetic relationship of individuals in the jth environment
## K[j,m] is the genetic relationship between individuals from the jth and the mth environment
## B[j] is the background genetic relationship of individuals in the jth environment
## B[j,m] is the background genetic relationship between individuals from the jth and the mth environment
## hg[j] and hb[j] are the genetic and background genetic heritability of a phenotypic trait in the jth environment respectively
## hg[j,m] and hb[j,m] are the genetic and background genetic heritability of a phenotypic trait between the jth and mth environment respectively
## tst[j].mat is testing index data for the jth environments

update.Sigma=function(theta,K1,K2,K3,K12,K21,K13,K31,K23,K32,B1,B2,B3,B12,B21,B13,B31,B23,B32,O1,O2,O3,O12,O21,O13,O31,O23,O32) {
	gtheta1=theta[1];  gtheta2=theta[2];  gtheta3=theta[3]; # var_g1/g2/g3
	gtheta12=theta[4]; gtheta13=theta[5]; gtheta23=theta[6]; # rho_g12/g13/g23
	btheta1=theta[7];   btheta2=theta[8];   btheta3=theta[9]; # var_b1/b2/b3
	btheta12=theta[10]; btheta13=theta[11]; btheta23=theta[12]; # rho_b12/b13/b23
	etheta1=theta[13];  etheta2=theta[14];  etheta3=theta[15]; # var_e1/e2/e3

  	gcov12=sqrt(gtheta1*gtheta2)*gtheta12; gcov13=sqrt(gtheta1*gtheta3)*gtheta13; gcov23=sqrt(gtheta2*gtheta3)*gtheta23;
	update.KG=rbind( cbind(gtheta1*K1, gcov12*K12, gcov13*K13),
                       cbind(gcov12*K21, gtheta2*K2, gcov23*K23),
			     cbind(gcov13*K31, gcov23*K32, gtheta3*K3 ) );

	bcov12=sqrt(btheta1*btheta2)*btheta12; bcov13=sqrt(btheta1*btheta3)*btheta13; bcov23=sqrt(btheta2*btheta3)*btheta23;
	update.KB=rbind( cbind(btheta1*B1, bcov12*B12, bcov13*B13),
                       cbind(bcov12*B21, btheta2*B2, bcov23*B23),
			     cbind(bcov13*B31, bcov23*B32, btheta3*B3 ) );

   	update.RE=rbind( cbind(etheta1*O1,O12,O13),
                       cbind(O21,etheta2*O2,O23),
			     cbind(O31,O32,etheta3*O3) );

  	new.Sigma=as.matrix(update.KG+update.KB+update.RE);
  	return(new.Sigma)
}

log_lik = function(Y,theta,Sigma){
	mu1=theta[16]; mu2=theta[17]; mu3=theta[18];
    	MU=c(rep(mu1,nrow(Y)),rep(mu2,nrow(Y)),rep(mu3,nrow(Y)));
  	log.LY = dmvn(as.vector(Y), MU, Sigma, log=T); # cal lnP
    	return(log.LY); 
}

is.between <- function(x, a, b) {
	x > a & x < b
}

log_prior = function(theta,Y,hg1_sq,hg2_sq,hg3_sq,hb1_sq,hb2_sq,hb3_sq,hg12,hg13,hg23,hb12,hb13,hb23){
  	gtheta1=theta[1];  gtheta2=theta[2];  gtheta3=theta[3]; # var_g1/g2/g3
	gtheta12=theta[4]; gtheta13=theta[5]; gtheta23=theta[6]; # rho_g12/g13/g23
	btheta1=theta[7];   btheta2=theta[8];   btheta3=theta[9]; # var_b1/b2/b3
	btheta12=theta[10]; btheta13=theta[11]; btheta23=theta[12]; # rho_b12/b13/b23
	etheta1=theta[13];  etheta2=theta[14];  etheta3=theta[15]; # var_e1/e2/e3
	mu.theta1=theta[16]; mu.theta2=theta[17]; mu.theta3=theta[18]; # mu1/mu2/mu3

  ## priors for G
  	df_g1=df_g2=df_g3=df_b1=df_b2=df_b3=df_e1=df_e2=df_e3=10; 
  	y1=Y[,1]; y2=Y[,2]; y3=Y[,3];
  	var_Y1=(sd(y1))^2; var_Y2=(sd(y2))^2; var_Y3=(sd(y3))^2;
  	S_g1=hg1_sq*var_Y1*(df_g1-2)/df_g1; 
  	S_g2=hg2_sq*var_Y2*(df_g2-2)/df_g2; 
	S_g3=hg3_sq*var_Y3*(df_g3-2)/df_g3;
  	prior.g1=dinvchisq(gtheta1,df=df_g1,scale=S_g1,log=T);
 	prior.g2=dinvchisq(gtheta2,df=df_g2,scale=S_g2,log=T);
   	prior.g3=dinvchisq(gtheta3,df=df_g3,scale=S_g3,log=T);

  	sd=0.2;
  	g12_mu=hg12*cor(y1,y2);   
 	prior.g12=dnorm(gtheta12,g12_mu,sd,log=T);
	g13_mu=hg13*cor(y1,y3);   
 	prior.g13=dnorm(gtheta13,g13_mu,sd,log=T);
	g23_mu=hg23*cor(y2,y3);   
 	prior.g23=dnorm(gtheta23,g23_mu,sd,log=T);
 
  ## priors for B  
  	S_b1=hb1_sq*var_Y1*(df_b1-2)/df_b1; 
  	S_b2=hb2_sq*var_Y2*(df_b2-2)/df_b2; 
	S_b3=hb3_sq*var_Y3*(df_b3-2)/df_b3;
  	prior.b1=dinvchisq(btheta1,df=df_b1,scale=S_b1,log=T); 
  	prior.b2=dinvchisq(btheta2,df=df_b2,scale=S_b2,log=T); 
	prior.b3=dinvchisq(btheta3,df=df_b3,scale=S_b3,log=T); 

  	b12_mu=hb12*cor(y1,y2);  
  	prior.b12=dnorm(btheta12,b12_mu,sd,log=T); 
	b13_mu=hb13*cor(y1,y3);  
  	prior.b13=dnorm(btheta13,b13_mu,sd,log=T); 
	b23_mu=hb23*cor(y2,y3);  
  	prior.b23=dnorm(btheta23,b23_mu,sd,log=T); 

  ## priors for E
  	S_e1=(1-hg1_sq-hb1_sq)*var_Y1*(df_e1-2)/df_e1; 
  	S_e2=(1-hg2_sq-hb2_sq)*var_Y2*(df_e2-2)/df_e2;
	S_e3=(1-hg3_sq-hb3_sq)*var_Y3*(df_e3-2)/df_e3; 
  	prior.e1=dinvchisq(etheta1,df=df_e1,scale=S_e1,log=T); 
  	prior.e2=dinvchisq(etheta2,df=df_e2,scale=S_e2,log=T); 
 	prior.e3=dinvchisq(etheta3,df=df_e3,scale=S_e3,log=T);

  ## priors for intercept 
	prior.mu1=dnorm(mu.theta1,0,100,log=T); prior.mu2=dnorm(mu.theta2,0,100,log=T); prior.mu3=dnorm(mu.theta3,0,100,log=T);

  	return( 	prior.g1+prior.g2+prior.g3+prior.g12+prior.g13+prior.g23+
          		prior.b1+prior.b2+prior.b3+prior.b12+prior.b13+prior.b23+
          		prior.e1+prior.e2+prior.e3+prior.mu1+prior.mu2+prior.mu3	 )
}

rw.prop = function(theta,Y,sd){
  	gtheta1=theta[1];  gtheta2=theta[2];  gtheta3=theta[3]; # var_g1/g2/g3
	gtheta12=theta[4]; gtheta13=theta[5]; gtheta23=theta[6]; # rho_g12/g13/g23
	btheta1=theta[7];   btheta2=theta[8];   btheta3=theta[9]; # var_b1/b2/b3
	btheta12=theta[10]; btheta13=theta[11]; btheta23=theta[12]; # rho_b12/b13/b23
	etheta1=theta[13];  etheta2=theta[14];  etheta3=theta[15]; # var_e1/e2/e3
	mu.theta1=theta[16]; mu.theta2=theta[17]; mu.theta3=theta[18]; # mu1/mu2/mu3

	y1=Y[,1]; y2=Y[,2]; y3=Y[,3];
  	var_Y1=(sd(y1))^2; var_Y2=(sd(y2))^2; var_Y3=(sd(y3))^2;

  repeat {
 	prop.g1=rnorm(1,gtheta1,sd); prop.b1=rnorm(1,btheta1,sd); prop.e1=var_Y1-prop.g1-prop.b1;      
	if ( prop.g1>0 & prop.b1>0 & prop.e1>0 ) break
  }

  repeat {
 	prop.g2=rnorm(1,gtheta2,sd); prop.b2=rnorm(1,btheta2,sd); prop.e2=var_Y2-prop.g2-prop.b2;      
	if ( prop.g2>0 & prop.b2>0 & prop.e2>0 ) break
  }

  repeat {
 	prop.g3=rnorm(1,gtheta3,sd); prop.b3=rnorm(1,btheta3,sd); prop.e3=var_Y3-prop.g3-prop.b3;      
	if ( prop.g3>0 & prop.b3>0 & prop.e3>0 ) break
  }

  repeat {
     prop.g12=rnorm(1,gtheta12,sd); prop.g13=rnorm(1,gtheta13,sd); prop.g23=rnorm(1,gtheta23,sd);
     if ( is.between(prop.g12,-1,1)==T & is.between(prop.g13,-1,1)==T & is.between(prop.g23,-1,1)==T) break
  }

    repeat {
     prop.b12=rnorm(1,btheta12,sd); prop.b13=rnorm(1,btheta13,sd); prop.b23=rnorm(1,btheta23,sd);
     if ( is.between(prop.b12,-1,1)==T & is.between(prop.b13,-1,1)==T & is.between(prop.b23,-1,1)==T) break
  }

	prop.mu1=rnorm(1,mu.theta1,sd); prop.mu2=rnorm(1,mu.theta2,sd); prop.mu3=rnorm(1,mu.theta3,sd);
  return( c( prop.g1,prop.g2,prop.g3,prop.g12,prop.g13,prop.g23,
		 prop.b1,prop.b2,prop.b3,prop.b12,prop.b13,prop.b23,
             prop.e1,prop.e2,prop.e3,prop.mu1,prop.mu2,prop.mu3) )
}

RWMH_chain = function(startvalue,iterations,Y,K1,K2,K3,K12,K21,K13,K31,K23,K32,B1,B2,B3,B12,B21,B13,B31,B23,B32,O1,O2,O3,O12,O21,O13,O31,O23,O32,
                      hg1_sq,hg2_sq,hg3_sq,hb1_sq,hb2_sq,hb3_sq,hg12,hg13,hg23,hb12,hb13,hb23) {
   chain = array(dim = c(iterations+1,length(startvalue)))
   chain[1,]=startvalue;
   for (i in 1:iterations){
    repeat {
    # do something
    curr.theta=chain[i,];
    prop.theta = rw.prop(curr.theta,Y,0.10);
    # exit if the condition is met
    Sigma.prop=update.Sigma(prop.theta,K1,K2,K3,K12,K21,K13,K31,K23,K32,B1,B2,B3,B12,B21,B13,B31,B23,B32,O1,O2,O3,O12,O21,O13,O31,O23,O32)
    if ( is.positive.definite(Sigma.prop)==T ) break
    }
    Sigma.now=update.Sigma(curr.theta,K1,K2,K3,K12,K21,K13,K31,K23,K32,B1,B2,B3,B12,B21,B13,B31,B23,B32,O1,O2,O3,O12,O21,O13,O31,O23,O32);
    #is.positive.definite(Sigma.now)==T;
    
    ratio=exp( log_lik(Y,prop.theta,Sigma.prop) + 
		   log_prior(prop.theta,Y,hg1_sq,hg2_sq,hg3_sq,hb1_sq,hb2_sq,hb3_sq,hg12,hg13,hg23,hb12,hb13,hb23) -
               log_lik(Y,curr.theta,Sigma.now) -
	         log_prior(curr.theta,Y,hg1_sq,hg2_sq,hg3_sq,hb1_sq,hb2_sq,hb3_sq,hg12,hg13,hg23,hb12,hb13,hb23) );
    probab = min(ratio,1);
    if (runif(1) < probab){ chain[i+1,] = prop.theta }
    else{chain[i+1,] = curr.theta}
  }
  return(chain) 
}


GBhat_3env=function(post.mean,y1,y2,y3,K1,K2,K3,K12,K21,K13,K31,K23,K32,B1,B2,B3,B12,B21,B13,B31,B23,B32,O1,O2,O3,O12,O21,O13,O31,O23,O32) {
	var_g1=post.mean[1]; var_g2=post.mean[2]; var_g3=post.mean[3]; 
	rho_g12=post.mean[4];rho_g13=post.mean[5];rho_g23=post.mean[6];  
	gcov12=sqrt(var_g1*var_g2)*rho_g12; gcov13=sqrt(var_g1*var_g3)*rho_g13; gcov23=sqrt(var_g2*var_g3)*rho_g23;
	KG=rbind( cbind(var_g1*K1, gcov12*K12, gcov13*K13),
          	    cbind(gcov12*K21, var_g2*K2, gcov23*K23),
	    	    cbind(gcov13*K31, gcov23*K32, var_g3*K3) );                         	
	var_b1=post.mean[7]; var_b2=post.mean[8]; var_b3=post.mean[9]; 
	rho_b12=post.mean[10];rho_b13=post.mean[11];rho_b23=post.mean[12];  
	bcov12=sqrt(var_b1*var_b2)*rho_b12; bcov13=sqrt(var_b1*var_b3)*rho_b13; bcov23=sqrt(var_b2*var_b3)*rho_b23;
	KB=rbind( cbind(var_b1*B1, bcov12*B12, bcov13*B13),
          		cbind(bcov12*B21, var_b2*B2, bcov23*B23),
	    		cbind(bcov13*B31, bcov23*B32, var_b3*B3) );                        
	var_e1=post.mean[13]; var_e2=post.mean[14]; var_e3=post.mean[15];
	RE=rbind( cbind(var_e1*O1,O12,O13),cbind(O21,var_e2*O2,O23),cbind(O31,O32,var_e3*O3) ); 

	varY=KG+KB+RE; 
	Ey1=rep(post.mean[16],length(y1)); Ey2=rep(post.mean[17],length(y2)); Ey3=rep(post.mean[18],length(y3));
	Y.diff=c( (y1.trn-Ey1), (y2.trn-Ey2), (y3.trn-Ey3) );
	Ghat=matrix(KG %*% solve(varY) %*% Y.diff,ncol=3); # size nx3
	Bhat=matrix(KB %*% solve(varY) %*% Y.diff,ncol=3); # size nx3
      res=cbind(Ghat,Bhat);
	return(res);
}

######################################## conduct cross-validated prediction #################################################################################
pred1=pred2=pred3=vector(mode="numeric",length=nrow(tst1.mat));
reps=50; 
y1hat=y2hat=y3hat=matrix(,ncol=reps,nrow=nrow(tst1.mat)); # save reps result

##setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores) #not to overload your computer
registerDoParallel(cl)

ptm <- proc.time() # star the clock
res=foreach( r=1:reps,.combine=cbind,.packages=c("matrixcalc","coda","geoR","mvnfast","foreach","doParallel","parallel","iterators")) %dopar% {
		tst1=tst1.mat[,r]; trn1=setdiff(seq(1,length(Y1)),tst1);
		tst2=tst2.mat[,r]; trn2=setdiff(seq(1,length(Y2)),tst2);  
		tst3=tst3.mat[,r]; trn3=setdiff(seq(1,length(Y3)),tst3); 
		y1.trn=Y1[trn1]; y2.trn=Y2[trn2]; y3.trn=Y3[trn3]; n.trn=length(trn1); 
		var1=round((sd(y1.trn))^2,4); var2=round((sd(y2.trn))^2,4); var3=round((sd(y3.trn))^2,4);				
		K1.trn=GK1[trn1,trn1]; K2.trn=GK2[trn2,trn2]; K3.trn=GK3[trn3,trn3]; 
		K12.trn=GK12[trn1,trn2]; K13.trn=GK13[trn1,trn3]; K23.trn=GK23[trn2,trn3];				
		B1.trn=B1[trn1,trn1]; B2.trn=B2[trn2,trn2]; B3.trn=B3[trn3,trn3];
		B12.trn=B12[trn1,trn2]; B13.trn=B13[trn1,trn3]; B23.trn=B23[trn2,trn3];
		O1.trn=O2.trn=O3.trn=diag(n.trn); O12.trn=O13.trn=O23.trn=matrix(0,nrow=n.trn,ncol=n.trn); 
	
		## get parameter estimation
		cor12=round(cor(y1.trn,y2.trn),2); cor13=round(cor(y1.trn,y3.trn),2); cor23=round(cor(y2.trn,y3.trn),2);

		start=c( 	hg1*var1,hg2*var2,hg3*var3,hg12*cor12,hg13*cor13,hg23*cor23,
	   			hb1*var1,hb2*var2,hb3*var3,hb12*cor12,hb13*cor13,hb23*cor23,
         			var1*(1-hg1-hb1),var2*(1-hg2-hb2),var3*(1-hg3-hb3),
        			mean(y1.trn),mean(y2.trn),mean(y3.trn) );
		nIter=100,000; burnIn=nIter*0.5;

		chain=RWMH_chain( startvalue=start, iterations=nIter, Y=cbind(y1.trn,y2.trn,y3.trn),
                  		K1=K1.trn, K2=K2.trn, K3=K3.trn, K12=K12.trn, K21=t(K12.trn),
					K13=K13.trn, K31=t(K13.trn), K23=K23.trn, K32=t(K23.trn), 
					B1=B1.trn, B2=B2.trn, B3=B3.trn, B12=B12.trn, B21=t(B12.trn), 
					B13=B13.trn, B31=t(B13.trn), B23=B23.trn, B32=t(B23.trn),
                  		O1=O1.trn, O2=O2.trn, O3=O3.trn, O12=O12.trn, O21=t(O12.trn), 
					O13=O13.trn, O31=t(O13.trn), O23=O23.trn, O32=t(O23.trn), 
                  		hg1_sq=hg1, hg2_sq=hg2, hg3_sq=hg3, hg12=hg12, hg13=hg13, hg23=hg23, 
					hb1_sq=hb1, hb2_sq=hb2, hb3_sq=hb3, hb12=hb12, hb13=hb13, hb23=hb23 );				

		chain_afterBurn=chain[-(1:(burnIn+1)),];
		mcmc_chain=mcmc(chain_afterBurn);
		post_mean=c();
		for (j in 1:ncol(mcmc_chain)) {
			post_mean[j]=mean(mcmc_chain[,j]); 
		}
		GB.hat=GBhat_3env( post_mean,y1.trn,y2.trn,y3.trn,K1.trn,K2.trn,K3.trn,K12.trn,t(K12.trn),K13.trn,t(K13.trn),K23.trn,t(K23.trn),
		       	 	B1.trn,B2.trn,B3.trn,B12.trn,t(B12.trn),B13.trn,t(B13.trn),B23.trn,t(B23.trn),
				 	O1.trn,O2.trn,O3.trn,O12.trn,t(O12.trn),O13.trn,t(O13.trn),O23.trn,t(O23.trn));
		Ghat=GB.hat[,1:3]; Bhat=GB.hat[,4:6]; # blup.g and blup.b from trn data
		K.trn=rbind( cbind(K1.trn,K12.trn,K13.trn),cbind(t(K12.trn),K2.trn,K23.trn),cbind(t(K13.trn),t(K23.trn),K3.trn) ); # dim(K.trn);
		B.trn=rbind( cbind(B1.trn,B12.trn,B13.trn),cbind(t(B12.trn),B2.trn,B23.trn),cbind(t(B13.trn),t(B23.trn),B3.trn) ); # dim(B.trn);

		## tst_trn data
		K1.tst.trn=GK1[tst1,trn1]; K2.tst.trn=GK2[tst2,trn2]; K3.tst.trn=GK3[tst3,trn3]; 
		K12.tst.trn=GK12[tst1,trn2]; GK21=t(GK12); K21.tst.trn=GK21[tst2,trn1];
		K13.tst.trn=GK13[tst1,trn3]; GK31=t(GK13); K31.tst.trn=GK31[tst3,trn1];
		K23.tst.trn=GK23[tst2,trn3]; GK32=t(GK23); K32.tst.trn=GK32[tst3,trn2];
		K.tst.trn=rbind( cbind(K1.tst.trn,K12.tst.trn,K13.tst.trn),
			     		cbind(K21.tst.trn,K2.tst.trn,K23.tst.trn),
			     		cbind(K31.tst.trn,K32.tst.trn,K3.tst.trn) ); # dim(K.tst.trn);
		Gtst.hat=matrix(K.tst.trn %*% solve(K.trn) %*% as.vector(Ghat),ncol=3); # dim(Gtst.hat)

		B1.tst.trn=B1[tst1,trn1]; B2.tst.trn=B2[tst2,trn2]; B3.tst.trn=B3[tst3,trn3];
		B12.tst.trn=B12[tst1,trn2]; B21=t(B12); B21.tst.trn=B21[tst2,trn1];
		B13.tst.trn=B13[tst1,trn3]; B31=t(B13); B31.tst.trn=B31[tst3,trn1];
		B23.tst.trn=B23[tst2,trn3]; B32=t(B23); B32.tst.trn=B32[tst3,trn2];
		B.tst.trn=rbind( cbind(B1.tst.trn,B12.tst.trn,B13.tst.trn),
			     		cbind(B21.tst.trn,B2.tst.trn,B23.tst.trn),
			     		cbind(B31.tst.trn,B32.tst.trn,B3.tst.trn) ); # dim(B.tst.trn);
		Btst.hat=matrix(B.tst.trn %*% solve(B.trn) %*% as.vector(Bhat),ncol=3); # dim(Btst.hat)
	
		Ey1=rep(post_mean[16],length(tst1)); Ey2=rep(post_mean[17],length(tst2)); Ey3=rep(post_mean[18],length(tst3));
		pred1=Ey1+Gtst.hat[,1]+Btst.hat[,1]; pred2=Ey2+Gtst.hat[,2]+Btst.hat[,2]; pred3=Ey3+Gtst.hat[,3]+Btst.hat[,3];
		return( cbind(pred1,pred2,pred3) );
	}
stopCluster(cl); #stop cluster
proc.time() - ptm   # stop the clock
	v=reps*3; # for 3 envs
	v1=c(seq(1,v,by=3)); v2=c(seq(2,v,by=3)); v3=c(seq(3,v,by=3));
	y1hat=res[,v1]; 
	y2hat=res[,v2];
	y3hat=res[,v3];
	
