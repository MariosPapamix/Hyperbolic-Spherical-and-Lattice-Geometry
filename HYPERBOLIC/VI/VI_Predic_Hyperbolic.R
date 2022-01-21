smpl_us_from_q <- function(tm, tsig, tus, tu_s, geom, bk_id = c(1,2)){
  # function to sample \alpha and z's from variational family
  #print(tus)
  ## sample \alpha
  a_smp = rnorm( 1, tm, tsig )
  
  ## sample latent coordinates
  N = dim( tus )[1]
  
  us_smp = array( NA, dim( tus ) )
  #print(us_smp)
  if (geom == "hyperbolic"){
    
    #us_smp[bk_id[1],] = tus[bk_id[1],]
    us_smp[bk_id[1],] = tus[bk_id[1],]
    #print(tus[i,])
    #print(tu_s[i])
    for (i in (1:N)[-bk_id[1]]){
      us_smp[i,] = sample_from_RHN( 1, tus[i,], tu_s[i])
    }
    us_smp[bk_id[2],] = tus[bk_id[2],]
    
  } else if (geom == "spherical"){
    
    us_smp[bk_id[1],] = tus[bk_id[1],]
    for (i in (1:N)[-bk_id[1]]){
      us_smp[i,] = rmovMF( 1, matrix( tus[i,] * tu_s[i], nrow=1) )
    }
    us_smp[bk_id[2],] = tus[bk_id[2],]
    
  } else if (geom == "discrete"){
    
    us_smp[bk_id[1],] = tus[bk_id[1],]
    for (i in (1:N)[-bk_id[1]]){
      #us_smp[i,] = rmovMF( 1, matrix( tus[i,] * tu_s[i], nrow=1) )
    }
    us_smp[bk_id[2],2] = tus[bk_id[2],2]
    
  } else { stop("Invalid choice for geom")}
  
  return( list( alpha = a_smp, us = us_smp ) )
}

calc_pij_vals <- function(alpha, us, geom){
  
  if ( (any( c("euclidean","hyperbolic", "spherical","discrete") == geom)) == FALSE ){
    stop("Invalid choice for geom")
  }
  
  N = dim(us)[1]
  ps = matrix(0, N, N)
  for (i in 1:(N-1)){
    for (j in (i+1):N){
      eta = alpha - sqrt((us[i,1]+ us[j,1])^2+(us[i,2]+ us[j,2])^2) 
      ps[i,j] = 1 / ( 1 + exp(-eta))
      ps[j,i] = ps[i,j]
    }
  }
  
  return(ps)
}



calc_post_pred <- function(tm, tsig, tus, ts_u, geom, M = 1000){
  # prm 1-4 are parameters of q
  # geom = "hyperbolic", "spherical" indicates the geometry
  # M is number of posterior samples
  
  ## create storage
  N = dim( tus )[1]
  ps_strg = array( NA, c(N,N,M) )
  
  for (m in 1:M){
    
    
    ## sample from variational family
    
    smp = smpl_us_from_q(tm[m], tsig[m], tus[,,m], ts_u[,m], geom)
    
    ## calculate connection probabilities from sample
    ps_strg[,,m] = calc_pij_vals(smp$alpha, smp$us, geom)
  }
  
  
  
  return( ps_strg )
}

preds<-calc_post_pred(tst$tm, tst$tsig, tst$tus, tst$ts_u, "hyperbolic", M = 1000)



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
  
  plot( y1dens, col='black', ylim=c(0, max(c(y0dens$y, y1dens$y))), xlim=c(0,1), xlab="Predictive probability of a link", main="VI-Stein Predict Fit" )
  lines( y0dens, col='red' ) 
  legend("topright", legend=c("Yij=0", "Yij=1"), col=c('red', 'black'), lwd=1)
}

dens_plot_pijs(Y,preds)
