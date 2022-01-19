# Brownian motion on Poincaré disk 
Random_Walk_Hyperbolic<-function(z,stepsize, numsteps){

#path<-c(z);
for (t in 1:numsteps) { 
  jitter=runif(1,0,1);
  dz=stepsize*complex(1,cos(2*pi*jitter),sin(2*pi*jitter));
  z<-(z+dz)/(Conj(dz)*z+1);
#  path<-c(path,z);
#  print(z)
  
}

#print(sqrt(Re(z)^2+Im(z)^2))

return(z)

}

# Brownian motion on the Sphere 


# Function based on movMF to simulate theta
Random_Walk_Spherical  <-function(Z,i){
  #n<-dim(Z)[1]

  Z_old<-Z
    if(i==2){
      theta<-0
    }else{
    theta<- (atan(Z[2]/Z[1])+rnorm(1,0,0.01))%%(2*pi)
    }
    phi<-(acos(Z[3])+rnorm(1,0,0.01))%%(2*pi)
  
  Z_old<-Z
  
 # print(theta)
#print(phi)
  
  if(i==2){
    Z_old[1]<-sin(phi)
    Z_old[2]<-0
    Z_old[3]<-cos(phi)
  }

    if(i>=3){
      Z_old[1]<-sin(phi)*cos(theta) 
      Z_old[2]<-sin(phi)*sin(theta) 
      Z_old[3]<-cos(phi)
    }
  
 #   print(Z_old)
    if(Z_old[1]>0 && i==2){
      Z<-Z_old
    }else if(Z_old[2]>0 && i==3){
      Z<-Z_old
    }else{
      Z<-Z_old
    }
    
    return(Z)
}



# Sphere to stereographic projection - Not needed 
  
sphere_to_stereographic<-function(output){
  
  output_stereographic <- matrix(1,P-1)
  
  output_stereographic[1,1]<-output[n+1,1]/(1-output[n+1,3]) 
  
  output_stereographic[1,2]<-output[n+1,2]/(1-output[n+1,3]) 
  
  return(output_stereographic)
}

# Brownian motion on Lattice 

Random_Walk_Discrete<-function(z,metric)
{
  xdir<-Re(z)
  ydir<-Im(z)
  

    r<-runif(1)
    if(r<=0.25) {xdir<-xdir+metric}
    if(r>0.25 && r<=0.5) {xdir<-xdir-metric}
    if(r>0.5 && r<=0.75) {ydir<-ydir +metric}
    if(r>0.75) {ydir<-ydir-metric}
    
  return(complex(real=xdir,imaginary=ydir))
}


Random_Walk_Continuous<-function(z,N){
  xdir<-Re(z)
  ydir<-Im(z)
  
  
   xdir<-xdir+rnorm(1,0,0.1)
   ydir<-ydir +rnorm(1,0,0.1)
  
  return(complex(real=xdir,imaginary=ydir))

 }

#print(path)

#plot(Re(path),Im(path),type="l",col="blue",asp=1,xlim=c(-1,1),ylim=c(-1,1));
#curve(sqrt(1-x^2),-1,1,col="red",add=TRUE);
#curve(-sqrt(1-x^2),-1,1,col="red",add=TRUE);