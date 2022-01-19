#ln distributions for the acceptance ratio 


lnprior_h<-function(Z,mu,sigma,n,i){
  
  prop_new<--h_distance(c(0,0),Z[i,])/sigma^2
  
  return(prop_new)
  
}


lnprior_s<-function(Z,mu,kappa,n,i){
  
  prop_new<-as.numeric(log(dvmf(Z[i], k=10, c(0,0,1))))
  
  return(prop_new)
  
}

lnprior_d<-function(zi,mu,sigma,n,k){

  prop_new<--pointDistance(mu, zi[k,],lonlat=FALSE)/(2*sigma^2)-log_S
  
  #print(prop_new)
  return(prop_new)
  
}