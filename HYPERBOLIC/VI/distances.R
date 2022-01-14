# Function for calculating distances in poincare disk, sphere and lattice.


#distance on Poincare disk

h_distance <- function(mu,zi) {
  #dist_h<-acosh(1+2*dist(rbind(mu,zi))^2/((1-dist(rbind(c(0,0),mu))^2)*(1-dist(rbind(c(0,0),zi))^2)))
#  print(dist_h)
  
  #print(mu)
  
  z2<-complex(real=Re(mu[1]),imaginary=Re(mu[2]))
  
  z1<-complex(real=Re(zi[1]),imaginary=Re(zi[2]))
  
  dist_h<-log(((abs(1-Conj(z1)*z2)+abs(z2-z1))/(abs(1-Conj(z1)*z2)-abs(z2-z1))))
  
  return(dist_h)
}

#Distance on the sphere


s_distance <- function(mu,zi) {
  
    dist_s<-acos( mu[1]*zi[1]+mu[2]*zi[2]+mu[3]*zi[3] )
  
  return(dist_s)
}


# Distance for the lattice |Re(x_i)-Re(x_j)|+|Im(x_i)-Im(x_j)|

l_distance <- function(mu,zi) {
  #dist_h<-acosh(1+2*dist(rbind(mu,zi))^2/((1-dist(rbind(c(0,0),mu))^2)*(1-dist(rbind(c(0,0),zi))^2)))
  #  print(dist_h)
  
  #print(mu)
  
#  dist_d<-abs(zi[1]-mu[1])+abs(zi[2]-mu[2])
  dist_d<-sqrt(abs(zi[1]-mu[1])^2+abs(zi[2]-mu[2])^2)

  return(dist_d)
}

