## script to calculate predictive distributions from the fitted models
## for each: input variational family parameters values and output a posterior summary
## to calculate predictive: sample \alpha and z's from variational family and calculate intersection probabilities pij's
source("sample_network_functions.R")

smpl_us_from_mcmc <- function(A,Z,geom){
  # function to sample \alpha and z's from variational family
  
  x<-sample.int(9000, 10000)
  
  ## sample \alpha
  a_smp = alpha[x]
  
  ## sample latent coordinates
  N = dim( tus )[1]
  us_smp = array( NA, dim( tus ) )
  if (geom == "hyperbolic"){
    
    for (i in (1:N)){
      us_smp[i,] = Z.post[[x]][i,]
    }
    
  } else if (geom == "spherical"){
    
    for (i in (1:N)){
      us_smp[i,] = Z.post[[x]][i,]
    }

  } else if (geom == "discrete"){
    
    for (i in (1:N)){
      #us_smp[i,] = rmovMF( 1, matrix( tus[i,] * tu_s[i], nrow=1) )
      us_smp[i,] =Z.post[[x]][i,]
    }
  } else { stop("Invalid choice for geom")}
  
  return( list( alpha = a_smp, us = us_smp ) )
}

calc_pij_vals <- function(alpha, us, geom){
  
  if ( (any( c("hyperbolic", "spherical","discrete") == geom)) == FALSE ){
    stop("Invalid choice for geom")
  }
  
  N = dim(us)[1]
  ps = matrix(0, N, N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      eta = alpha - dist_geom( us[i,], us[j,], geom ) 
      ps[i,j] = 1 / ( 1 + exp(-eta))
      ps[j,i] = ps[i,j]
    }
  }
  
  return(ps)
}

## get ys

calc_post_pred <- function(A,Z, geom, M = 1000){
  # prm 1-4 are parameters of q
  # geom = "hyperbolic", "spherical" indicates the geometry
  # M is number of posterior samples
  
  ## create storage
  N = dim( tus )[1]
  ps_strg = array( NA, c(N,N,M) )
  
  for (m in 1:M){
    
    ## sample from variational family
    smp = smpl_us_from_mcmc(A,Z, geom)
    
    ## calculate connection probabilities from sample
    ps_strg[,,m] = calc_pij_vals(smp$alpha, smp$us, geom)
  }
  return( ps_strg )
}

## get ys

calc_preds_mcmc <- function(utrc, atrc, rng, geom){
  nsmp = length(rng)
  N = dim(utrc)[1]
  preds = array( NA, c(N,N,nsmp) )
  for (a in 1:nsmp){
    preds[,,a] = calc_pij_vals(atrc[rng[a]], utrc[,,rng[a]], geom)
  }
  return( preds )
}

get_upper_tris <- function(ys, preds){
  id_uppertri = upper.tri( ys )
  id_y1 = id_uppertri & (ys==1) # T if upper triangle and yij=1
  id_y0 = id_uppertri & (ys==0) # T if upper triangle and yij=0
  
  ## extract upper triangles for each
  M = dim(preds)[3]
  p_y1 = matrix(NA, sum(id_y1), M)
  p_y0 = matrix(NA, sum(id_y0), M)
  for (m in 1:M){
    ptmp = preds[,,m]
    p_y1[,m] = ptmp[id_y1]
    p_y0[,m] = ptmp[id_y0]
  }
  return( list(y1s=p_y1, y0s=p_y0) )
}

dens_plot_pijs <- function(ys, preds){
  preds = get_upper_tris( ys, preds )
  
  y1dens = density( preds$y1s )
  y0dens = density( preds$y0s )
  
  plot( y1dens, col='black', ylim=c(0, max(c(y0dens$y, y1dens$y))), xlim=c(0,1), xlab="Predictive probability of a link", main="" )
  lines( y0dens, col='red' ) 
  legend("topright", legend=c("Yij=0", "Yij=1"), col=c('red', 'black'), lwd=1)
}

plot_pij_preds <- function( ys, preds, mean=FALSE, ylims=NA ){
  
  id_uppertri = upper.tri( ys )
  id_y1 = id_uppertri & (ys==1) # T if upper triangle and yij=1
  id_y0 = id_uppertri & (ys==0) # T if upper triangle and yij=0
  
  ## extract upper triangles for each
  M = dim(preds)[3]
  p_y1 = matrix(NA, sum(id_y1), M)
  p_y0 = matrix(NA, sum(id_y0), M)
  for (m in 1:M){
    ptmp = preds[,,m]
    p_y1[,m] = ptmp[id_y1]
    p_y0[,m] = ptmp[id_y0]
  }
  
  if (mean==FALSE){
    ## plot on histogram
    xmin = min( preds )
    xmax = max( preds )
    col1 = rgb(1,0,0,0.5)
    col2 = rgb(0,0,1,0.5)
    par(mfrow=c(1,1))
    if (is.na(ylims)==TRUE){
      hist( p_y0, xlim=c(xmin, xmax), col=col1, breaks=30, main="", xlab="probability" )
    } else {
      hist( p_y0, xlim=c(xmin, xmax), ylim=ylims, col=col1, breaks=30, main="",xlab="probability")
    }
    hist( p_y1, add=T, col=col2, breaks=30)
    legend( "topright", legend = c("y=0", "y=1"), col = c(col1, col2), pch=19)
  } else if (mean==TRUE){
    p_y0 = rowMeans(p_y0)
    p_y1 = rowMeans(p_y1)
    
    ## plot on histogram
    xmin = min( c(p_y0, p_y1) )
    xmax = max( c(p_y0, p_y1) )
    col1 = rgb(1,0,0,0.5)
    col2 = rgb(0,0,1,0.5)
    par(mfrow=c(1,1))
    if (is.na(ylims)==TRUE){
      hist( p_y0, xlim=c(xmin, xmax), col=col1, breaks=30, main="", xlab="probability" )
    } else {
      hist( p_y0, xlim=c(xmin, xmax), ylim=ylims, col=col1, breaks=30, main="",xlab="probability")
    }
    hist( p_y1, add=T, col=col2, breaks=30)
    legend( "topright", legend = c("y=0", "y=1"), col = c(col1, col2), pch=19)
  }
}
