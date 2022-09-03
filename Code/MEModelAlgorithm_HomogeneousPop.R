library(MTM)

## The code below is used to fit multi-environment (ME) model for same individuals across different environments.
## The ME model is fitted by using R package MTM MTM (de los Campos and Grüneberg, 2016) with 20,000 iterations and the first 10,000 samples as burn-in. 
## y is phenotype data of multiple environments
## WK is the genetic relationship of individuals estimated by weighted kernel
## B is the background genetic relationship of individuals that is not explained by genetic markers

fm <- MTM( Y = y,
           K = list( list( K = WK, COV = list( type = 'UN', df0=ncol(y), S0=diag(ncol(y)) ) ), 
                     list( K = B,  COV = list( type = 'UN', df0=ncol(y), S0=diag(ncol(y)) ) )),
                     resCov = list( type = 'DIAG', S0 = rep(1,ncol(y)), df0 = rep(1,ncol(y))),
                     nIter = 20000, burnIn = 10000, thin = 1, saveAt = 'model_fitting');
