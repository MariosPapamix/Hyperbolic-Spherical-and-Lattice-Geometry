## R code to calculate frechet means using python
library('reticulate') # for calling python in R
library('spacyr')
#py_install("geomstats", pip = TRUE)
#use_python("/usr/bin/python") # path to python (may be different on your machines)
#spacy_initialize(model = "en_core_web_lg")
#reticulate::use_condaenv("spacy_condaenv", required = TRUE)
use_python("C:/Users/SuperM/AppData/Local/r-miniconda/envs/r-reticulate/python.exe")
#fuse_python("C:/Users/SuperM/anaconda3/python.exe")
source_python('frechet_mean.py') # contains function to calculate the mean

## function is calc_frechet_mean(dat, d=2, geom, maxit), where dat is (Nxdim) matrix, d is dimension of manifold (=2 in our examples), geom = "sphere" or "poincare_disk", maxit is maximum iterations for algorithm (I set to 500)

## EXAMPLE in the sphere
geom = "sphere"
maxit = 500
dat = matrix( c(0,0,1, 0, 1/sqrt(2), 1/sqrt(2)), ncol=3, byrow=TRUE )
fm_sphere = calc_frechet_mean(dat, 2, geom, maxit) 


##MCMC sphere

x_sphere_mcmc<-Z.post[[1]][9000:10000,]
y_sphere_mcmc<-Z.post[[2]][9000:10000,]
z_sphere_mcmc<-Z.post[[3]][9000:10000,]
geom = "sphere"
maxit = 5000
fm_sphere_mcmc=matrix(0,15,3)
for(i in 1:15){
  dat=cbind(x_sphere_mcmc[,i],y_sphere_mcmc[,i],z_sphere_mcmc[,i])
  fm_sphere_mcmc[i,] = calc_frechet_mean(dat, 2, geom, maxit) 
}

library("rgl")
spheres3d(0,0,0,lit=FALSE,color="white")
spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")

spheres3d(fm_sphere_mcmc[,1],fm_sphere_mcmc[,2],fm_sphere_mcmc[,3],col=c(1:15),radius=0.1)


##VI sphere

geom = "sphere"
maxit = 5000
fm_sphere_vi=matrix(0,15,3)
for(i in 1:15){
  dat=cbind(tst$tu[i,1,501:1000],tst$tu[i,2,501:1000],tst$tu[i,3,501:1000])
  fm_sphere_vi[i,] = calc_frechet_mean(dat, 2, geom, maxit) 
}

library("rgl")
spheres3d(0,0,0,lit=FALSE,color="white")
spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")

spheres3d(fm_sphere_vi[,1],fm_sphere_vi[,2],fm_sphere_vi[,3],col="red",radius=0.01)



#fm_poincare_vi[,2]<--fm_poincare_vi[,2]


## EXAMPLE in the poincare disk
#geom = "poincare_disk"
#maxit = 5000
#dat = matrix( c(0,.1, .2, .4), ncol=2 )
#fm_poincare = calc_frechet_mean(dat, 2, geom, maxit) 


##MCMC poincare

x_hyperbolic_mcmc<-Z.post[[1]][9000:10000,]
y_hyperbolic_mcmc<-Z.post[[2]][9000:10000,]
geom = "poincare_disk"
maxit = 5000
fm_poincare_mcmc=matrix(0,34,2)
for(i in 1:34){
  dat=cbind(x_hyperbolic_mcmc[,i],y_hyperbolic_mcmc[,i])
  fm_poincare_mcmc[i,] = calc_frechet_mean(dat, 2, geom, maxit) 
}



##VI poincare

geom = "poincare_disk"
maxit = 5000
fm_poincare_vi=matrix(0,34,2)
for(i in 1:34){
  dat=cbind(tst$tus[i,1,1001:2000],tst$tus[i,2,1001:2000])
  fm_poincare_vi[i,] = calc_frechet_mean(dat, 2, geom, maxit) 
}

fm_poincare_vi[,2]<--fm_poincare_vi[,2]


par(mfcol = c(1, 2))
plot(fm_poincare_mcmc,col=c(1:34),pch = 16,xlab="x-axis",ylab="y-axis",main="HM-MCMC centroids of latent positions")
text(fm_poincare_mcmc,col="black")
plot(fm_poincare_vi,col=c(1:34),pch = 16,xlab="x-axis",ylab="y-axis",main="VI centroids of latent positions")
text(fm_poincare_vi,col="black")


