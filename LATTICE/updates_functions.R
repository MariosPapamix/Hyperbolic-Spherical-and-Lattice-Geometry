#Resolve identifiability issues in Spherical case

rotmat <- function( thtx, thty, thtz ){
  matrix( c( cos(thtz)*cos(thty), sin(thtz)*cos(thty), - sin(thty), cos(thtz)*sin(thty)*sin(thtx) - sin(thtz)*cos(thtx), sin(thtz)*sin(thty)*sin(thtx) + cos(thtz)*cos(thtx), cos(thty)*sin(thtx), cos(thtz)*sin(thty)*cos(thtx) + sin(thtz)*sin(thtx), sin(thtz)*sin(thty)*cos(thtx) - cos(thtz)*sin(thtx), cos(thty)*cos(thtx)), ncol=3, byrow=FALSE )
}

ident_sphere<- function(us, bk_id = c(1,2)){

    us_bk1 = us[bk_id[1],]
    us_bk2 = us[bk_id[2],]
    dist_bk = acos( sum( us_bk1 * us_bk2 ) )
    
    ## bk1 mapped to c(0,0,1), bk2 mapped to c(a,0,b) (calc a and b)
    b = cos( dist_bk )
    a = sqrt( 1 - b^2 )
    
    ## calculate parameters of transformation
    thtx = atan( us_bk1[2] / us_bk1[3] )
    thty = atan( - us_bk1[1] / (us_bk1[2]*sin(thtx) + us_bk1[3]*cos(thtx)) )
    thtz = atan( (us_bk2[3]*sin(thtx) - us_bk2[2]*cos(thtx)) / (us_bk2[1]*cos(thty) + sin(thty)*(us_bk2[2]*sin(thtx) + us_bk2[3]*cos(thtx) ) ) ) 
    
    ## apply transform to coordinates    
    R = rotmat( thtx, thty, thtz )
    us_t = t( R %*% t(us) )
    
    if(us_t[2,1]<0){
      
      for(i in 1:dim(us_t)[1]){
        
        us_t[i,1]<--us_t[i,1] 
        
      }
      
    }
    
    if(us_t[3,2]<0){
      
      for(i in 1:dim(us_t)[1]){
        
        us_t[i,2]<--us_t[i,2] 
        
      }
            
    }
    
    
    
    return( list(Z = us_t) )

}

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
  
  #b<-0.5*(1-Conj(z_c[1])*z_c[2])/(z_c[2]-z_c[1])
  
  #if(Re(b*(z_c[2]-a)/(1-Conj(a)*z_c[2]))<0){
  #  b<--b
  #}
  
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
  
  #reflection<-Z[3,2]
  
  #if(reflection<0){
  #    Z[3:34,2]<--Z[3:34,2]
    
#  }
  
#  print(Z)
  
  return(Z)
  
}

#Resolve identifiability issues in Lattice case **********

ident_lat<-function(Z,n){
  
  x<-Z[1,]
  
  for( i in 1:n) { Z[i,]<-Z[i,]-x } #translation
  
  x1<-Z[2,1]
  x2<-Z[2,2]
  x3<-Z[1,]
  x4<-Z[2,]

  
  x<-rbind(Z[1,],Z[2,])
  
  for(i in 1:n){  
    x5<-Z[i,1]
    Z[i,1]<-x5*dist(x)/x2+Z[i,2]*dist(x)/x1 
    Z[i,2]<-x5*dist(x)/x1-Z[i,2]*dist(x)/x2
  } #rotation
  
  for( i  in 2:n){
    if(Z[2,1]<0){
      Z[i,]<--Z[i,]
    }
  }#reflection
  
  for( i  in 3:n){
    if(Z[3,2]<0){
      Z[i,2]<--Z[i,2]
    }
  }#reflection
  
  metric<-Z[2,1]-Z[1,1]
  
  #print(Z)
  #print(metric)
  
  return(list(Z=Z,metric=metric))
  
##}

}
#Update latent position each time

Z.up<-function(Y,Z,alpha,n,metric){
  
  x=2
  #mu.z=c(0,0)
  mu.z=c(0,0,0)
  sd.z=10
  ##update Z
  Z_old<-Z
  
  #print(n)
  
  for(i in x:n){
  #stepsize<-sample(c(0.05), size=1, replace=TRUE, prob=c(1))
 # print(stepsize)
  #temp<-Random_Walk_Hyperbolic(complex(real=Z[i,1],imaginary=Z[i,2]),stepsize,1)
  temp<-Random_Walk_Spherical(Z[i,],i)
  #temp<-Random_Walk_Discrete(Z[i,1],imaginary=Z[i,2]),metric)
  #temp<-Random_Walk_Continuous(complex(real=Z[i,1],imaginary=Z[i,2]),1)
  
  
  ##For hyperbolic, continuous and lattice
  
  #Z_old[i,]<-Z[i,]
  ##if(i==2){
  ##  Z[i,1]<-abs(Re(temp))
  ##  Z[i,2]<-0
  ##}
  ##else if(i==3 && Z[3,2]<0){
  ##  Z[i,1]=Re(temp)
  ##  Z[i,2]=abs(Im(temp))
  ##}else{
  ##Z[i,1]=Re(temp)
  ##Z[i,2]=Im(temp)
  ##}
  #random walk Z
  }

  #print(Z)
  #print(alpha)


    for(k in 1:(n-1)){
      lnew<-0
      loldz<-0
      for(j in (k+1):n){
        
        #eta<-alpha-h_distance(Z[k,],Z_old[j,])
        eta<-alpha-s_distance(Z[k,],Z_old[j,])
        #eta<-alpha-l_distance(Z[k,],Z_old[j,])
        loglik1<-eta*Y[k,j]-log(1+exp(eta))
        
        #eta<-alpha-h_distance(Z_old[k,],Z_old[j,])
        eta<-alpha-s_distance(Z_old[k,],Z_old[j,])
        #eta<-alpha-l_distance(Z_old[k,],Z_old[j,])
        loglik2<-eta*Y[k,j]-log(1+exp(eta))
        
        lnew<-lnew+loglik1
        loldz<-loldz+loglik2
      }
      
      #print(Z)
      #print(alpha)
      
      #print(Z)
      #print(Z_old)

     #hr<-lnew-loldz+lnprior_d(Z,mu.z,sd.z,n,k)-lnprior_d(Z_old,mu.z,sd.z,n,k)
      hr<-lnew-loldz+lnprior_s(Z,mu.z,sd.z,n,k)-lnprior_s(Z_old,mu.z,sd.z,n,k)
     # hr<-lnew-loldz+log(dnorm(Z[k,],c(0,0),10))-log(dnorm(Z_old[k,],c(0,0),10))
      
     #print(lnew)
     #print(loldz)
     #print(exp(hr))
    
    if( runif(1)>exp(hr)) {
      Z[k,]<-Z_old[k,]
    }

  }
  
  lnew<-0
  
  for(k in 1:(n-1)){
    for(j in (k+1):n){
        
        #eta<-alpha-h_distance(Z_old[k,],Z_old[j,])
        eta<-alpha-s_distance(Z_old[k,],Z_old[j,])
        #eta<-alpha-l_distance(Z_old[k,],Z_old[j,])
        loglik1<-eta*Y[k,j]-log(1+exp(eta))
        lnew<-lnew+loglik1

      
    }
  }
  
  list(Z1=Z,lik=lnew)
  
}

########################################

#Update alpha each time

alpha.up<-function(Y,Z,alpha,adelta,n,a.a=1,a.b=1,lold){
  ##update alpha
  stepsize<-sample(c(0.05), size=1, replace=TRUE, prob=c(1))
  alphanew<-alpha+rnorm(1, 0, stepsize) 
  
  #print(n)
  
  log_lik<-0
  
  for(k in 1:(n-1)){
    for(j in (k+1):n){

        
        #eta<-alphanew-h_distance(Z[k,],Z[j,])
        eta<-alphanew-s_distance(Z[k,],Z[j,])
        #eta<-alphanew-l_distance(Z[k,],Z[j,])
        loglik1<-eta*Y[k,j]-log(1+exp(eta))
        log_lik<-log_lik+loglik1
      
    }
  }
  #print(log_lik)
  
  lnew<-log_lik
  
  hr<-exp( lnew-lold )*( dnorm(alphanew,0,10)/dnorm(alpha,0,10) ) 
 
  #print(hr)
  
  if(runif(1)<hr){

    alpha<-alphanew
    lold<-lnew
    
  }
  
  list(alpha=alpha,lik=lold)  
  
}

########################################
