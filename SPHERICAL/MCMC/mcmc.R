##set random seed
#set.seed(1)

#install.packages("igraph")
#install.packages("smacof") # for MDS in sphere
#install.packages("pracma") # needed for erf()
#install.packages("mvtnorm")
#install.packages("plotrix")
#install.packages("rgl")
#install.packages("plot3D")
#install.packages("movMF")
#install.packages("Rfast")
#install.packages("raster")
#install.packages("igraphdata")
#install.packages("mva")    
#install.packages("netrankr")

library(Directional)
library(Rfast)
library(igraph)
library(smacof) # for MDS in sphere
library(pracma) # needed for erf()
library(mvtnorm)
library(plotrix)
library(rgl)
library(plot3D)
library(movMF) # needs to be one or the other - could be causing issues
library(raster)
library(netrankr)
library(network) 
library(ergm)
library(igraphdata)
library(latentnet)


source("initialization.R") 
source("distances.R") 
source("Random_Walks.R")
source("updates_functions.R")
source("sample_network_functions.R")
source("lnprior.R")



##library mva provides cdmscale, used to get starting values
#library(mva)    

#print(1111)

##read Karate Club data

#Y<-karate$adjacency

##read florentine 

#$(florentine_m) 
florentine_m <- delete_vertices(florentine_m,which(degree(florentine_m)==0))

Y<-as.matrix(as_adjacency_matrix(florentine_m))

##read Monks

#data(sampson)
#Y<-as.matrix(samplike)

##take out unconnected actors
zro<-(1:dim(Y)[1])[ Y%*%rep(1,dim(Y)[1])==0 ] 
if(length(zro)>0) {Y<-Y[-zro,-zro]}


n<-dim(Y)[1]    #number of nodes
k<-2            #dimension of latent space


## initial values of Z and alpha

D<-dist(Y)

Z<-matrix(0,nrow=n,ncol=k)

#Z<-init_us_poincare(Y) 
Z<-init_us_sphere(Y) 
#Z<-init_us_lattice(Y) 

metric<-0

#Z<-ident_hyp(Z,n)
Z<-ident_sphere(Z) 

Z<-as.matrix(Z[[1]])
#x<-ident_lat(Z,n)
#Z<-x$Z
#metric<-x$metric
#metric<-metric/100

#for(i in 1:n){
#  Z[i,1]<-Z[i,1]-(Z[i,1]%%metric)
#  Z[i,2]<-Z[i,2]-(Z[i,2]%%metric)
#}

#log_S<-0

#for(i in 1:100){
#  for(j in 1:100){
#    log_S<-log_S+2*exp(-pointDistance(c(0,0), c(i*metric,j*metric),lonlat=FALSE)/(2*10)^2) 
#  }
#}

#log_S<-log(1+log_S)

#alpha<-init_alpha_poincare(Y,Z)
alpha<-init_alpha_sphere(Y,Z)
#alpha<-init_alpha_lattice(Y,Z)

print(alpha)

loglik<-0

#  print(Z)

#initial likelihood

for(k in 1:(n-1)){
  for(j in (k+1):n){
    if(k!=j){
      eta<-alpha-s_distance(Z[k,],Z[j,])
      loglik<-loglik+eta*Y[k,j]-log(1+exp(eta))
    }
  }
}

k<-2

adelta<-.5          #params for
zdelta<-.2          #proposal distribution

nscan<-10^3         #number of scans   
odens<-10^3         #save output every odens step

aca<-acz<-0          #keep track of acceptance rates
Alpha<-alpha        #keep track of alpha


#lik<-lpY(Y,Z,alpha,-100000,n) #keep track of alpha and likelihood 

Z.post<-list()    #keep track of positions
for(i in 1:k){ Z.post[[i]]<-t(Z[,i]) }

## MCMC

lik<-loglik

lik_mat<-matrix(0,1000,1)

for(ns in 1:1000){
 
  print(ns)
  print(lik)

  tmp<-Z.up(Y,Z,alpha,n,metric)                     #update z's

  acz<-acz+1/odens

  Z<-tmp$Z1
  
  lik<-tmp$lik

  tmp<-alpha.up(Y,Z,alpha,adelta,n,a.a=2,a.b=1,lik)  #update alpha
  
  if(tmp$alpha!=alpha) { 
    aca<-aca+1/odens
    alpha<-tmp$alpha 
    lik<-tmp$lik
  }
  
  Alpha<-c(Alpha,alpha)
  acz<-aca<-0 
  for(i in 1:2){ Z.post[[i]]<-rbind(Z.post[[i]],t(Z[,i])) } 
  
  lik_mat[ns]<-lik
  
  ##plot(Z.post[[1]][i,],Z.post[[2]][i,])
  
}

Post<-list(Z=Z.post,Alpha=Alpha,Lik=lik)
dput(Post,"ff.m.post") #output to a file

##plot results, for k=2
Zp<-Post$Z
#Lik<-Post$Lik
Alpha<-Post$Alpha


##plot likelihood and alpha
#par(mfrow=c(1,2))
#plot(Lik,type="l")
#plot(Alpha,type="l")

##plot positions and confidence sets
#par(mfrow=c(1,1))
#par(mar=c(3,3,1,1))
#par(mgp=c(2,1,0))


