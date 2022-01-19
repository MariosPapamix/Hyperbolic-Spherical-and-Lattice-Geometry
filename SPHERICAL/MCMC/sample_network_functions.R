library(pracma) # needed for erf()
library(mvtnorm)
library(plotrix)
library(rgl)
library(plot3D)
library(movMF) # needs to be one or the other - could be causing issues

plot_graph <- function(us, ys, geom="euclidean"){

    N = dim(us)[1] # number of nodes

    if (geom == "euclidean"){
        minx = min(us[,1])
        maxx = max(us[,1])
        miny = min(us[,2])
        maxy = max(us[,2])
        
        plot( us[,1], us[,2], pch=19, xlim=c(minx, maxx), ylim=c(miny, maxy) )
        text( us[,1], us[,2], labels=1:N, pos=1)
        for (i in 1:(N-1)){
            for (j in (i+1):N){
                if (ys[i,j]==1){
                    segments( us[i,1], us[i,2], us[j,1], us[j,2])
                }
            }
        }
    } else if (geom == "hyperbolic"){       
        plot( us[,1], us[,2], pch=19, xlim=c(-1,1), ylim=c(-1,1) )
        text( us[,1], us[,2], labels=1:N, pos=1)
        for (i in 1:(N-1)){
            for (j in (i+1):N){
                if (ys[i,j]==1){
                    segments( us[i,1], us[i,2], us[j,1], us[j,2])
                    ## NOTE: plotting doesn't reflect geometry (lines are Euclidean)
                }
            }
        }
    }
    
    ## else if (geom == "spherical"){

    ##     spheres3d(0,0,0,lit=FALSE,color="white")
    ##     spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")
    ##     spheres3d(netsmp$us[,1], netsmp$us[,2], netsmp$us[,3], col='black',radius=.02)
    ##     text3d(netsmp$us[,1], netsmp$us[,2], netsmp$us[,3], 1:N, pos=3)

    ##     for (i in 1:(N-1)){
    ##         for (j in (i+1):N){
    ##             if (ys[i,j]==1){
    ##                 segments3d( x0=us[i,1], y0=us[i,2], z0=us[i,3], x1=us[j,1], y1=us[j,2], z1=us[j,3], add=TRUE, col='black', lwd=2)
    ##             }
    ##         }
    ##     }
    ## }
    
}

plot_coords_sphere <- function(us){
    N = dim(us)[1]
    spheres3d(0,0,0,lit=FALSE,color="white")
    spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")
    spheres3d(netsmp$us[,1], netsmp$us[,2], netsmp$us[,3], col='black',radius=.02)
    text3d(us[,1], us[,2], us[,3], 1:N, pos=3)
}



################################### SPHERICAL #########################################

calc_gamma_1 <- function(theta){
    c(cos(theta[2]), sin(theta[2])*cos(theta[1]), sin(theta[2])*sin(theta[1]) )
}
calc_gamma_2 <- function(theta){
    c(- cos(theta[3])*sin(theta[2]), cos(theta[3])*cos(theta[2])*cos(theta[1]) - sin(theta[3])*sin(theta[1]), cos(theta[3])*cos(theta[2])*sin(theta[1]) + sin(theta[3])*cos(theta[1]))
}
calc_gamma_3 <- function(theta){
    c(sin(theta[3])*sin(theta[2]), - sin(theta[3])*cos(theta[2])*cos(theta[1]) - cos(theta[3])*sin(theta[1]), - sin(theta[3])*cos(theta[2])*sin(theta[1]) + cos(theta[3])*cos(theta[1]) )
}

calc_tht_from_gamma <- function(g1, g2, g3){
    tht2 = acos( g1[1] )
    tht1 = acos( g1[2] / sin(tht2) )
    tht3 = asin( g3[1] / sin(tht2) )
    return( c(tht1, tht2, tht3) )
}

spherical_dist <- function( x, y ){ # distance on a sphere
    acos( sum(x*y) )
}

#################################### LATTICE ##########################################

## To do...

#######################################################################################
#################################### Shell function ###################################
#######################################################################################

smpl_us <- function( N, d, uprms, geom ){
    if (geom == "hyperbolic" ){
        ##us = sample_from_wrapped_norm(N, d, uprms$sig, uprms$mn)$z_exp # points on the hyperbolic
        ##yspoinc = map_hyperboloid_to_poincare(ushyp) # maps points to poincare (better for plotting)
        us = sample_from_RHN( N, uprms$mu, uprms$sig)
    } else if (geom == "spherical"){
        ## us = rkent( N, uprms$kappa, uprms$mn, uprms$beta )
        us = rmovMF( N, matrix( uprms$mn * uprms$k, nrow=1) )
    } else if (geom == "euclidean"){
        us = rmvnorm( N, mean = uprms$mu, sigma = uprms$Sigma )
    }
    return( us )
}

dist_geom <- function( x, y, geom ){
    if (geom == "hyperbolic"){
        ## dst = lorentz_dist(x,y)
        dst = hypdistpoinc(x,y)
    } else if (geom == "spherical"){
        dst = spherical_dist( x, y )
    } else if (geom == "euclidean"){
        dst = sqrt( sum( (x - y)^2 ) )
    }
    return( dst )
}

sample_network <- function(N, d, alpha, uprms, geom, tau=NA){ 
                                        # N is number of nodes
                                        # d is dimension of latent space
                                        # alpha is base rate for edges
                                        # geom = "hyperbolic", "spherical", "lattice" or "euclidean"
                                        # uprms is list of distn prms for Us
                                        # tau is threshold for connection in RGG graph (leave as NA to
                                        # sample from the model)

    ## create storage
    ys = matrix(0, N, N) # stores edges
    ps = matrix(0, N, N) # stores probabilities
    ds = matrix(0, N, N) # stores distances
    
    ## sample latent coordinates
    us = smpl_us( N, d, uprms, geom )    
    
    ## sample edges
    for (i in 1:(N-1)){
        for (j in (i+1):N){
            ds[i,j] = dist_geom( us[i,], us[j,], geom )
            if (is.na(tau) == FALSE ){
                if (ds[i,j] < tau ){ # have a thresholded version to make testing easier
                    ps[i,j] = 1
                    ys[i,j] = 1
                } else{
                    ps[i,j] = 0
                    ys[i,j] = 0
                }
            } else {
                eta = alpha - ds[i,j]
                ps[i,j] = 1 / (1 + exp(-eta))
                ys[i,j] =  rbinom(1,1, ps[i,j] )
            }
            ## populate summetric
            ds[j,i] = ds[i,j]
            ps[j,i] = ps[i,j]
            ys[j,i] = ys[i,j]
        }
    }

    return( list(ys=ys, ps=ps, us=us, ds=ds) )
}

