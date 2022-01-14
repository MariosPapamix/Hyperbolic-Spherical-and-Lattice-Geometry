#ln distributions for the acceptance ratio 


lnprior_h<-function(Z,mu,sigma,n,i){
  
  prop_new<--h_distance(c(0,0),Z[i,])^2/(2*sigma^2)
  
  return(prop_new)
  
}