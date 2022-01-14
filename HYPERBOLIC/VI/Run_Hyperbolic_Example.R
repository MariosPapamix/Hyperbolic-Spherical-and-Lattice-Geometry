source('bbvi_hyperbolic.R')
###################### check code ##################################
set.seed(12375)

library(hydra)

N = 34
d = 2
alpha = 1.5
geom = "hyperbolic"
uprms = list(sig = .5, mu=c(0,0)) 
netsmp = sample_network(N, d, alpha, uprms, geom, tau=1.2)

par(mfrow=c(1,1))


source("initialisation.R") 
source("sample_network_functions.R")



##library mva provides cdmscale, used to get starting values
#library(mva)    
library(igraphdata)


##read Karate Club data

Y<-karate$adjacency

##read florentine 

#data(florentine_m) 
#florentine_m <- delete_vertices(florentine_m,which(degree(florentine_m)==0))

#Y<-as_adjacency_matrix(florentine_m)

##read Monks

#samplike

##take out unconnected actors
#zro<-(1:dim(Y)[1])[ Y%*%rep(1,dim(Y)[1])==0 ] 
#if(length(zro)>0) {Y<-Y[-zro,-zro]}


n<-dim(Y)[1]    #number of nodes
k<-2   #dimension of latent space

N<-n

## initial values of Z and alpha

D<-dist(Y)

Z<-matrix(0,nrow=n,ncol=k)


Z<-init_us_poincare(Y)

Z<-ident_hyp(Z,n)

#print(Z)
print(n)
#Z<-ident_sphere(Y) 
#x<-ident_lat(Y)
#Z<-x$Z
#metric<-x$metric
#metric<-metrix/100

alpha<-init_alpha_poincare(Y,Z)
#alpha<-init_alpha_sphere(Y,Z)
#alpha<-init_alpha_lattice(Y,Z)

print(alpha)
print(1111)

loglik<-0

for(k in 1:(n-1)){
  for(j in (k+1):n){
    eta<-alpha-h_distance(Z[k,],Z[j,])
    loglik<-loglik+eta*Y[k,j]-log(1+exp(eta))
  }
}
print(loglik)
plot_graph(Z, Y, "hyperbolic")

print(Z)
print(alpha)
print(uprms$sig)
print(uprms$mu)

## ## parameters
prms_p = list( us=Z, m=alpha, sig=uprms$sig, mu=uprms$mu, su = uprms$sig )

## ## run BBVI 
S = 5
maxrep = 2000
st_tim <- Sys.time()
Rprof()
tst = hyperbolic_BBVI(Y, prms_p, S, maxrep)
Rprof(NULL)
end_tim <- Sys.time()

summ = summaryRprof()

#pal = colorRampPalette(colors = c('grey', 'red'))
#cols = pal(maxrep+1)
#col_nodes = rainbow(N)

#plt_rng = 1:(maxrep+1)
#plot_output(tst, plt_rng, cols, col_nodes)
