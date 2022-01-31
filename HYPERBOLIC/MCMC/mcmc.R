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
#install.packages("remotes")
#remotes::install_github("ahoundetoungan/PartialNetwork")
#install.packages("hydra")


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
#library(netrankr)
library(network) 
library(ergm)
library(igraphdata)
#library(latentnet)
library(circular)

source("initialization.R") 
source("distances.R") 
source("Random_Walks.R")
source("updates_functions.R")
#source("sample_network_functions.R")
source("lnprior.R")


##read Karate Club data

Y<-as.matrix(karate$adjacency)

#Y<-as.matrix(read.table("ff.m.dat"))

#Y<-Y[order(colSums(Y),decreasing = TRUE),order(colSums(Y),decreasing = TRUE)]

#x<-Y[c(15,12,3,4,5,6,7,8,9,10,11,1,13,14,2),c(15,12,3,4,5,6,7,8,9,10,11,1,13,14,2)]

#Y<-x

#Y<-matrix(c(0,1,1,1,0,0,1,0,0),3,3)

#g <- sample_pa(100)

#Y<- as.matrix(as_adjacency_matrix(g))

##take out unconnected actors
zro<-(1:dim(Y)[1])[ Y%*%rep(1,dim(Y)[1])==0 ] 
if(length(zro)>0) {Y<-Y[-zro,-zro]}

aca<-acz<-0          #keep track of acceptance rates

n<-dim(Y)[1]    #number of nodes
k<-2            #dimension of latent space


## initial values of Z and alpha

##read Karate Club data

#Y<-karate$adjacency

##read florentine 

#data(florentine_m) 
#Y <- delete_vertices(florentine_m,which(degree(florentine_m)==0))

#Y<-as_adjacency_matrix(Y)
##read Monks

#samplike

##take out unconnected actors



n<-dim(Y)[1]    #number of nodes
k<-2   #dimension of latent space

N<-n

## initial values of Z and alpha

D<-dist(Y)

Z<-matrix(0,nrow=n,ncol=k)


Z<-init_us_poincare(Y)

Z<-ident_hyp(Z,n)

#print(Z)
#print(n)
#Z<-ident_sphere(Y) 
#x<-ident_lat(Y)
#Z<-x$Z
#metric<-x$metric
#metric<-metrix/100

alpha<-init_alpha_poincare(Y,Z)

loglik<-0


for(k in 1:(n-1)){
  for(j in (k+1):n){
      eta<-alpha-h_distance(Z[k,],Z[j,])
      loglik<-loglik+eta*Y[k,j]-log(1+exp(eta))
  }
}
print(loglik)

#Sys.sleep(100000)

k<-2

adelta<-.5          #params for
zdelta<-.2          #proposal distribution

odens<-10         #save output every odens step

aca<-acz<-0          #keep track of acceptance rates
Alpha<-alpha        #keep track of alpha


Z.post<-list()    #keep track of positions
for(i in 1:k){ Z.post[[i]]<-t(Z[,i]) }

## MCMC

lik<-loglik

nscan<-10^4        #number of scans 
lik_mat<-matrix(0,nscan,1)

start_time = Sys.time()



for(ns in 1:nscan){
 
  print(ns)
  
  tmp<-Z.up(Y,Z,alpha,n,metric,lik)                     #update z's

  if(tmp$Z1!=Z){
  #acz<-acz+1
  Z<-tmp$Z1
  lik<-tmp$lik
  }
  
  
  
  

  tmp<-alpha.up(Y,Z,alpha,adelta,n,a.a=2,a.b=1,lik)  #update alpha
  
  if(tmp$alpha!=alpha) { 
    aca<-aca+1
    alpha<-tmp$alpha 
    lik<-tmp$lik
  }
  

  Alpha<-c(Alpha,alpha)
  
  for(i in 1:2){ Z.post[[i]]<-rbind(Z.post[[i]],t(Z[,i])) } 
  
  lik_mat[ns]<-lik
  
  print(lik)
  print(alpha)
  
}

acz<-acz/nscan

aca<-aca/nscan

Post<-list(Z=Z.post,Alpha=Alpha,Lik=lik)
dput(Post,"ff.m.post") #output to a file

##plot results, for k=2
Zp<-Post$Z
#Lik<-Post$Lik
Alpha<-Post$Alpha

end_time = Sys.time()

time<-end_time - start_time

##plot likelihood and alpha
#par(mfrow=c(1,2))
#plot(Lik,type="l")
#plot(Alpha,type="l")

##plot positions and confidence sets
#par(mfrow=c(1,1))
#par(mar=c(3,3,1,1))
#par(mgp=c(2,1,0))


#acf(Alpha)
