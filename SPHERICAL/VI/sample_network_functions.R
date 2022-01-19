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

#######################################################################################
################################### HYPERBOLIC ########################################
## have functions for wrapped (Nagano et al) and Riemmanian/max entropy Gaussian

hypdistpoinc <- function(x,y){

    acosh( 1 + ( 2 * sum((x-y)^2) ) / ( (1-sum(x^2))*(1-sum(y^2))  ) )
}

lorentz_product <- function(x, y){
    len = length(x)
    lprod = sum( x[2:len]*y[2:len] ) - x[1]*y[1]
    return( lprod )
}

lorentz_dist <- function( x, y ){
    acosh( - lorentz_product(x,y) )
}

PT <- function(mu0, mu, v){
    alpha = - lorentz_product(mu0, mu)
    pt = v + ( lorentz_product(mu - alpha*mu0, v) / (alpha+1) ) * (mu0 + mu)
    return(pt)
}

invPT <- function(mu, mu0, u){
    inv_pt = PT(mu, mu0, u)
    return(inv_pt)
}

lorentz_norm <- function( x ){
    return( sqrt( lorentz_product(x, x) ) )
}

exp_map <- function(x, mu){
    xnorm = lorentz_norm(x)
    expm = cosh( xnorm )*mu + (sinh( xnorm )/xnorm)*x 
    return( expm )
}

inv_exp_map <- function(x, mu){
    alpha = - lorentz_product(mu,x)
    inv_expm = (acosh( alpha) / sqrt( alpha^2 - 1 )) * (x - alpha*mu)
    return( inv_expm )
}

map_hyperboloid_to_poincare <- function(z){
    len = dim(z)[2]
    return( z[,2:len]/(1+z[,1]) )
}

map_poincare_to_hyperboloid <- function(z){
    len = dim(z)[2] # number of dimensions
    denom = 1 - rowSums( z^2 )
    z_hyp = cbind( (1 + rowSums( z^2 ))/denom, 2*z/denom )
    return( z_hyp )
}

sample_from_wrapped_norm <- function(n, mu, sig){
    ## set mu0 (always the same)
    d = length(mu) - 1
    mu0 = c(1, rep(0, d) )

    ## sample normal
    vtilde = matrix(rmvnorm( n, c(0,0), sig ), nrow=n)
    v = cbind( rep(0,n), vtilde )

    ## apply PT
    u = matrix(NA, ncol=d+1, nrow=n)
    for (i in 1:n){ u[i,] = PT(mu0, mu, v[i,]) }

    ## 4) apply exponential map
    z = matrix(NA, ncol=d+1, nrow=n)
    for (i in 1:n){ z[i,] = exp_map(u[i,], mu) }

    return(list(v=v, u_pt=u, z_exp=z))
}

pdf_wrappednorm <- function(z, mu, cov, d=2, log=TRUE){
    ## note: not vectorised, this is Nagano distribution
    mu0 = c(1, rep(0, d) )

    ## apply inverse maps
    u = inv_exp_map(z, mu)
    v = invPT(mu, mu0, u)
    unorm = lorentz_norm( u )
    
    ## calculate density
    pdf = dmvnorm( v[-1], c(0,0), cov ) #note: 1st entry is deterministic = 0
    pdf = pdf * ( unorm / sinh(unorm) )^(d-1)
    if (log==TRUE){ pdf = log(pdf) }

    return(pdf)
}

add_mobius <- function(x, y){
    num = (1 + 2*sum(x*y) + sum(y^2) )*x + (1 - sum(x^2))*y
    denom = 1 + 2*sum(x*y) + sum(x^2)*sum(y^2)
    return( num / denom )
}

lam_z <- function( z ){
    2 / (1 - sum( z^2 ) )
}

exp_map_add_mob <- function(x, mu){
    ## exponential map expressed as a function of mobius addition (may be equivalent to above)
    lam_mu = lam_z( mu )
    exp_mu = add_mobius( mu, tanh( lam_mu * sqrt( sum( x^2 ) ) / 2 ) * x / sqrt( sum( x^2 ) ) )
    return( exp_mu )
}

z_hypgauss <- function(sig){
    2 * pi * sqrt( pi/2 ) * sig * exp( (sig^2)/2 ) * erf( sig / sqrt(2) )
}

pdf_hypgauss <- function(x, mu, sig){
    exp( - hypdistpoinc(x, mu)^2 / (2 * sig^2) ) / z_hypgauss(sig)
}

rho_r <- function(r, sig){
    exp( - r^2 / (2 * sig^2) ) * sinh(r) / z_hypgauss(sig)
}

sample_from_RHN <- function( N , mu, sig){
    ## RHN = Riemann hyperbolic Normal
    ## use algorithm described in Mathieu et al (2019)
    ## take Gamma rejection sampler (only looking at d=2 so dimensionality isn't an issue)
    ## NOTE: may return NaN AR is mu is taken 'too close' to the boundary
    
    ## calculate the envelope (fixed throughout)
    M = gamma(2) * (sig^2) * exp( ( sig + 1)^2 / 2 ) / z_hypgauss(sig) 
    

    smpl = matrix(0, N, 2) # initialise storage
    ind = 1 # index for storing samples

    ## run sampler
    while (ind <= N){
        
        ## sample uniform (for accept reject)
        u = runif(1)

        ## sample from proposal - use wrapped normal with same mu
        angl = runif(1, 0, 2*pi) # sample a uniform direction
        dir = sqrt( runif(1) ) # sample a magnitude
        alpha = dir * c( cos(angl), sin(angl) ) # uniform point in 2d disk      
        r = rgamma(1, shape=2, scale=sig) # sample a gamma magnitude
        
        ## convert to coordinates (ie not polar )
        x_poinc = exp_map_add_mob( r * alpha / lam_z(mu), mu )
               
        ## calculate ratio of prop to distn - this is wrong!!
        ar = rho_r(r, sig) / dgamma(r, 2, sig) 

        ## accept with prob
        
        if(is.na(ar/M)==FALSE){

        if ( runif(1) < (ar / M) ){            
            smpl[ind,] = x_poinc # store sample (in Poincare disk)
            ind = ind + 1 # update indicator for samples
        }
        
        }

    }

    return( smpl )
}

## PREVIOUS VERSION: keep just in case
## sample_from_hypgauss <- function( N = 10, mu, sig, mHat = 2, sig_prop=NA ){
##     ## mHat is the upper bound for the envelope (is updated as part of the sampler)
##     ## known as 'empirical supremum rejection sampling'
##     ## note: I tried a uniform sampler and couldn't get it to work properly...
   
##     smpl = matrix(0, N, 2) # initialise vector of
##     ind = 1 # index for storing samples
##     if( is.na(sig_prop)==TRUE ){ sig_prop = 2*sig*diag(2) } # inflate target sigma for proposal

##     mu_hp = map_poincare_to_hyperboloid( matrix(mu, nrow=1) ) # mean on hyperboloid
    
##     ## run sampler
##     while (ind <= N){
        
##         ## sample uniform
##         u = runif(1)

##         ## sample from proposal - use wrapped normal with same mu
##         x = sample_from_wrapped_norm( 1, mu_hp, sig_prop )$z_exp
##         x_poinc = map_hyperboloid_to_poincare(x) # map to poincare
##         ## calculate ratio of prop to distn
##         ar = pdf_hypgauss(x_poinc, mu, sig) / pdf_wrappednorm( x, mu_hp, sig_prop, log=FALSE )
        
##         ## accept with prob
##         if ( runif(1) < (ar / mHat) ){
##             smpl[ind,] = x_poinc # store sample (in Poincare disk)
##             ind = ind + 1 # update indicator for samples

##             ##if ( ind %% 10 == 0 ){ print(ind) }
##         }

##         ## update envelope
##         mHat = max( mHat, ar )
##     }

##     return( list(smpl = smpl, mHat = mHat) )
## }


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

lattice_dist <- function( x, y ){ # distance on a sphere
    return(abs(x[1]-y[1])+abs(x[2]-y[2]))
}


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
      else if (geom == "euclidean"){
        dst = abs((x[1]-y[1])+abs(x[2]-y[2]))
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

