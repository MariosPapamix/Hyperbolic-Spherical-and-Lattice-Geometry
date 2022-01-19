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
  #print(us)
  #print(us[1,])
  N = dim(us)[1]
  #print(N)
  ps = matrix(0, N, N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      #print(us[j,])
      eta = alpha - s_distance( us[i,], us[j,] ) 
      ps[i,j] = 1 / ( 1 + exp(-eta))
      ps[j,i] = ps[i,j]
    }
  }
  
  return(ps)
}

## get ys

calc_post_pred <- function(A,Z, geom, M = 9000){
  # prm 1-4 are parameters of q
  # geom = "hyperbolic", "spherical" indicates the geometry
  # M is number of posterior samples
  
  ## create storage
  N = length( Z$Z[[1]][1,])
  #print(N)
  ps_strg = array( NA, c(N,N,1000) )
  
  for (m in 1:1000){
    
    ### sample from variational family
    #smp = smpl_us_from_mcmc(A,Z, geom)
    
    ## calculate connection probabilities from sample
    #ps_strg[,,m] = calc_pij_vals(smp$alpha, smp$us, geom)
    pr_Z<-cbind(Z$Z[[1]][m+M,],Z$Z[[2]][m+M,])
    #print(pr_Z)
    pr_A<-A[m+M]
    #print(pr_A)
    ps_strg[,,m] = calc_pij_vals(pr_A,pr_Z, geom)
  }
  return( ps_strg )
}

preds<-calc_post_pred(Alpha,Post,"spherical")

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


dens_plot_pijs(Y,preds)