#Resolve identifiability issues in Hyperbolic case

ident_hyp<-function(Z,n){
  
  
  z_c<-matrix(0,nrow=n,ncol=1)
  
  z_c_new<-matrix(0,nrow=n,ncol=1)
  
  a<-complex(real = Z[1,1], imaginary = Z[1,2])
  
  z_c[1]<-a
  
  z_c[2]<-complex(real = Z[2,1], imaginary = Z[2,2])
  
  b<-sqrt((cosh(h_distance(Z[1,],Z[2,]))-1)/
            (cosh(h_distance(Z[1,],Z[2,]))+1))*
    (Conj(z_c[1])*z_c[2]-1)/(z_c[2]-z_c[1])
  
  
  
  for(i in 1:n){
    z_c[i]<-complex(real = Z[i,1], imaginary = Z[i,2])
    z_c_new[i]<-b*(z_c[i]-a)/(1-Conj(a)*z_c[i])
    Z[i,1]<-Re(z_c_new[i])
    Z[i,2]<-Im(z_c_new[i])
  }
  
  if(Z[2,1]<0){
    for(i in 1:n){
      Z[i,1]<--Re(z_c_new[i])
    }
  }
  
  if(Z[3,1]<0){
    for(i in 1:n){
      Z[i,2]<--Im(z_c_new[i])
    }
  }
  
  
  return(Z)
  
}


#Update latent position each time

Z.up<-function(Y,Z,alpha,n,metric,lold){
  
  mu.z=c(0,0)
  sd.z=0.1
  
  Z_old<-Z
  
  
  
  for(i in 2:n){
    stepsize<-sample(c(0.01), size=1, replace=TRUE, prob=c(1))
    
    temp<-Random_Walk_Hyperbolic(complex(real=Z[i,1],imaginary=Z[i,2]),stepsize,1)
    
    ##For hyperbolic
    if(h_distance(c(0,0),c(Re(temp),Re(Im(temp))))<10){
      if(i==2){
        Z[i,1]<-abs(Re(temp))
      }
      else if(i==3 && Z[3,2]<0){
        Z[i,1]=Re(temp)
        Z[i,2]=abs(Im(temp))
      }else{
        Z[i,1]=Re(temp)
        Z[i,2]=Im(temp)
      }
    }
  }
  #random walk Z
  
  for(k in 1:(n-1)){
    lnew<-0
    loldz<-0
    for(j in (k+1):n){
      if(k!=j){
        eta<-alpha-h_distance(Z[k,],Z_old[j,])
        loglik1<-eta*Y[k,j]-log(1+exp(eta))
        
        eta<-alpha-h_distance(Z_old[k,],Z_old[j,])
        loglik2<-eta*Y[k,j]-log(1+exp(eta))
        
        lnew<-lnew+loglik1
        loldz<-loldz+loglik2
      }
    }
    
    hr<-lnew-loldz+lnprior_h(Z,mu.z,0.1,n,k)-lnprior_h(Z_old,mu.z,0.1,n,k)
    
    if( runif(1)<exp(hr)) {
      Z_old[k,]<-Z[k,]
    }
    
  }
  
  lnew<-0
  
  for(k in 1:(n-1)){
    for(j in (k+1):n){
      eta<-alpha-h_distance(Z[k,],Z[j,])
      loglik1<-eta*Y[k,j]-log(1+exp(eta))
      lnew<-lnew+loglik1
    }
  }
  
  #  x<-0
  
  #  for(i in 1:n){
  #    x<-x+lnprior_h(Z,mu.z,sd.z,n,i)-lnprior_h(Z_old,mu.z,sd.z,n,i)
  #  }
  #  hr<-lnew-lold+x
  
  #  if( runif(1)<exp(hr)) {
  #      Z_old<-Z
  #  }
  
  return(list(Z1=Z_old,lik=lnew))
  
}

########################################

#Update alpha each time

alpha.up<-function(Y,Z,alpha,adelta,n,a.a=1,a.b=1,lold){
  ##update alpha
  #stepsize<-sample(c(0.05), size=1, replace=TRUE, prob=c(1))
  alphanew<-alpha+runif(1,-0.5,0.5) 
  
  log_lik<-0
  
  for(k in 1:(n-1)){
    for(j in (k+1):n){
      eta<-alphanew-h_distance(Z[k,],Z[j,])
      loglik1<-eta*Y[k,j]-log(1+exp(eta))
      log_lik<-log_lik+loglik1
    }
  }
  
  lnew<-log_lik
  
  hr<-exp( lnew-lold )*( dnorm(alphanew,0,1)/dnorm(alpha,0,1) ) 
  
  if(runif(1)<hr){
    
    alpha<-alphanew
    lold<-lnew
    
  }
  
  list(alpha=alpha,lik=lold)  
  
}

########################################
