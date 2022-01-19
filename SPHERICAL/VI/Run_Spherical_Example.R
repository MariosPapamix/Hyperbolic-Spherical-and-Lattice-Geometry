source('bbvi_spherical_vmf.R')
###################### run code and plot output  ##################################

library(Directional)
library(netrankr)

 set.seed(1210)
 N = 15
 d = 3
 alpha = 1 # keep this here
 uprms = list(k=4, mn=c(0,1/sqrt(2),1/sqrt(2))) ## sample from von mises fisher
 netsmp = sample_network(N, d, alpha, uprms, "spherical", tau=.6) 
 
 hist(netsmp$ds)

 clear3d() # just to visualise the network
 plot_coords_sphere(netsmp$us)
 g <- graph_from_adjacency_matrix( netsmp$ys, mode="undirected" )
 par(mfrow=c(1,1))
 plot( g, vertex.labels=1:N )

 
## ## populate prm_p
# prms_p = list( us=netsmp$us, m=alpha, sig=2, k=uprms$k, mu=uprms$mn )
 
 data(florentine_m) 
 florentine_m <- delete_vertices(florentine_m,which(degree(florentine_m)==0))
 
 Y<-as_adjacency_matrix(florentine_m)
 
 prms_p = list( us=netsmp$us, m=alpha, sig=2, k=2, mu=1 )

## ## ok - the problem is that the coordinates aren't guaranteed to be sampled at the center - could calculate this from the initialisation??
## ## bit hacky but it would solve this issue!

## ## run BBVI
 S = 10
 maxrep = 1000
 st_tim <- Sys.time()
 Rprof()
 tst = spherical_BBVI(Y, prms_p, S, maxrep)
 prms_p = tst$prms_p
 tst = tst$hst
 Rprof(NULL)
 end_tim <- Sys.time()

 summ = summaryRprof()

## ## plot the output 
 pal = colorRampPalette(colors = c('grey', 'red'))
 cols = pal(maxrep+1)
 col_nodes = rainbow(N)

 plt_rng = 1:(maxrep+1)
 plot_traces(tst, plt_rng, cols, col_nodes)
 plot_spherical_coords(tst$tu, plt_rng, col_nodes)
 spheres3d(prms_p$mu[1], prms_p$mu[2], prms_p$mu[3], col='red', radius=0.03) 

## ## ##########################################
## ## sanity check isometry:
 us = (rmovMF( 10, 5*c(0,0,-1), alpha = 1))
 us_t = apply_isometry(us, t(c(0,0,1)))$us_t

## ## need to close the plot! that's why it was looking awful
 spheres3d(0,0,0,lit=FALSE,color="white")
 spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")
 spheres3d(us[,1], us[,2], us[,3], col='black', radius=0.03) 
 text3d(us[,1], us[,2], us[,3], 1:10, pos=3)
 spheres3d(us_t[,1], us_t[,2], us_t[,3], col='red', radius=0.03) 
 text3d(us_t[,1], us_t[,2], us_t[,3], 1:10, pos=3, col='red')

 ds = matrix(NA, 10, 10)
 ds_t = matrix(NA, 10, 10)
 for (i in 1:10){
     for (j in 1:10){
         if (i != j ){
         ds[i,j] = acos( sum( us[i,]*us[j,] ) )
         ds_t[i,j] = acos( sum( us_t[i,]*us_t[j,] ) )
         }
     }
 }
## round(ds - ds_t, 2) ## should be 0's
