## ok, this is code to estimate parameters of hyperbolic latent space network model
##rm(list= ls()) # removes all from workspace

source("sample_network_functions.R") # function to sample a network
library(igraph)
library(hydra)

## the model is p(y|u,alpha) p(alpha|m,sig) p(u| mu, s) 
## note: coordinates assumed to follow a Riemannian hyperbolic normal
## throughout, rhn is Riemann hyperbolic normal
## the variational family is q(alpha|tm,tsig) q(u|tus, ts)

##################### init storage #############################

init_strg <- function(N, maxItr){    
    hst = list(
        ## variational parameters
        tm = rep(NA, maxItr+1),
        tsig = rep(NA, maxItr+1),
        tus = array(NA, c(N, 2, maxItr+1)), # mean for coords
        tus_plr = array(NA, c(N, 2, maxItr+1)), # mean for coords
        ts_u = matrix(NA, N, maxItr+1), # variance for coords
        ELBO = rep(NA, maxItr+1) # store ELBO at each iteration
    )
    return( hst )
}

##################### extra useful functions #############################

samp_us_from_rhn <- function(S, tus, ts_u, bk_id){

    dim_tus = dim(tus)
    us_smp = array(0, c(dim_tus, S))
    ids_nonbk = (1:dim_tus[1])[-bk_id[c(1,3)]]
    for (i in ids_nonbk){       
        us_smp[i,,] = t( sample_from_RHN( S, tus[i,], ts_u[i]) )
    }
    # sample 3rd bookstein coordinate so 2nd entry positive
    #idkp = FALSE
    #while( idkp == FALSE ){
      #  us_smp[bk_id[3],,] = t( sample_from_RHN( S, tus[bk_id[3],], ts_u[bk_id[3]]) )
     #   if ( all(us_smp[bk_id[3],2,] > 0) ){ idkp = TRUE }
   # }
    #print(1111)
    us_smp[bk_id[3],2,]=abs(us_smp[bk_id[3],2,])    
    # constrain 1st and 2nd bk coordinates 
    us_smp[bk_id[1],,] = tus[bk_id[1],] 
    us_smp[bk_id[2],2,] = tus[bk_id[2],2] 

    return( us_smp )
}

samp_us_from_rhn_ind <- function(S, ind, tus, ts_u){

    N = dim(tus)[1]
    d = dim(tus)[2]
    
    us_smp = array(0, c(N, d, S))
    ## sample for index ind
    us_smp[ind,,] = t( sample_from_RHN( S, tus[ind,], ts_u[ind]) )
    ## add in other coordinates (note fixed at current us value)
    us_smp[-ind,,] = array( tus[-ind,], c(N-1, d, S))
    
    return( us_smp )

}

samp_from_q_ind <- function(S, tmcur, tsigcur, tuscur, tu_sigcur, ind, us){
    ## this function only samples new coordinate for node ind
    alph_s = rnorm( S, mean=tmcur, sd=tsigcur )     
    us_s = samp_us_from_rhn_ind(S, ind, tuscur, tu_sigcur)
    return( list(alph_s=alph_s, us_s=us_s) )
}

samp_from_q <- function(S, tmcur, tsigcur, tuscur, tu_sigcur, bk_id ){
    alph_s = rnorm( S, mean=tmcur, sd=tsigcur )
    us_s = samp_us_from_rhn(S, tuscur, tu_sigcur, bk_id)
    return( list(alph_s=alph_s, us_s=us_s) )
}

us_to_plr <- function(us){

    if( is.vector(us) == TRUE ){ us = matrix( us, nrow=1) }
    
    ## radii
    rs = sqrt( rowSums(us^2) )

    ## angl phi
    us1pos = us[,1] > 0
    us2pos = us[,2] > 0
    ##need to add on multiple of pi so that phi \in [0, 2pi) (multiple depends on xy quadrant)
    id_u10 = (us[,1] == 0)
    
    phi = rep(0, dim(us)[1] )
    phi = atan( us[,2]/us[,1] ) + pi*( (!us1pos & us2pos) | (!us1pos & !us2pos) ) + (2*pi)*( us1pos & !us2pos )

    phi[id_u10] = 0 # define angle = 0 when 1st coord = 0

    return( cbind(rs, phi) )
}

plr_to_us <- function(plrs){

    if( is.vector(plrs) == TRUE ){ plrs = matrix( plrs, nrow=1) }
    
    mat = cbind( cos(plrs[,2]), sin(plrs[,2]) )
    us = sweep( mat, 1, plrs[,1], '*') # multiple by the radii
    
    return(us)
}

## functions to convert between parameters and unconstrained updates
phi_to_phistr <- function( x ){ log( x ) - log( 2*pi - x ) }
phistr_to_phi <- function( x ){ 2*pi / (1 + exp( -x )) }
r_to_rstr <- function( x ){ log( x ) - log( 1 - x ) }
rstr_to_r <- function( x ){ 1 / (1 + exp( -x )) }

pi_to_r <- function( x ){ log( x ) - log( pi - x ) }
r_to_pi <- function( x ){ pi / (1 + exp( -x )) }

######################## log densities ##################################

log_py <- function(obs, us, alpha){ 
    
    ## note: should be easy to vectorise this

    ll = 0
    N = dim(us)[1]
    
    

    for (i in 1:(N-1)){
        for (j in (i+1):N){
            eta = alpha - hypdistpoinc( us[i,], us[j,] )
            
            ll = ll + obs[i,j]*eta - log(1 + exp(eta))
        }
    }
    
    
    #print(ll)

    return( ll )
}

log_py_ui <- function(obs, us, alpha, ind){ ## note: should be easy to vectorise
    ## return only the terms in log p(y | u, alpha ) which contain node ind
    
    ll = 0
    N = dim(us)[1]

    for (j in (1:N)[-ind]){ # loop over all node pairs
        eta = alpha - hypdistpoinc( us[ind,], us[j,] )
        ll = ll + obs[ind,j]*eta - log(1 + exp(eta))
    }
    return( ll )
}

log_rhn_pdf_ui <- function(u, mu, s){
    ## log p(u_i | mu, s)
    lpdf = - log( 2*pi*sqrt(pi/2) ) - log(s) - (s^2 /s) - log(erf(s/sqrt(2))) - ( hypdistpoinc( u, mu )^2 / (2 * s^2) )
    
    #print(555)
    #print(lpdf)
    #print(mu)

    #print(666)
    #return( lpdf )
    return( lpdf )
}

log_post_alpha <- function(alpha, m, sig){
    ## log p(alpha | m, sig )
    lp = dnorm(alpha, mean=m, sd=sig, log=TRUE )
    return( lp )
}

log_post_ui <- function(obs, us, alpha, mu, s_u, ind){
    ## calculate log p_i( y | u, alpha), log p(u_i | mu, s_u )
   lp = log_py_ui(obs, us, alpha, ind) + log_rhn_pdf_ui(us[ind,], mu, s_u )
    return( lp )
}

log_post_us <- function(us, mu, s_u){
    ## calculate sum_i log p(u_i | mu, s_u )
    N = dim(us)[1]   
    lp = 0
    #print(9999)
    #print(mu)
    #print(10101010)
    for ( i in (1:N) ){ lp = lp + log_rhn_pdf_ui(us[i,], mu, s_u) }
    return( lp )
}

log_q_alpha <- function(alpha, tm, tsig){
    ## calculate log q(alpha | tm, tsig )
    lq = dnorm(alpha, mean=tm, sd=tsig, log=TRUE )
    return( lq )
}

log_q_ui <- function(u, tu, ts_u){
    ## calculate log q( ui | tui, tsi )
  #print(2312)
  #print(tu)
  #print(3212)
    lq = log_rhn_pdf_ui(u, tu, ts_u)
    return( lq )
}

log_q_us <- function(us, tus, ts_u, bk_id){
    ## calculate sum log q( ui | tui, tsi )
    N = dim(us)[1]
    lq = 0

    for ( i in (1:N)[-bk_id][1] ){
     # print(2222)
     # print(tus[i,])
     # print(8888)
        lq = lq + log_rhn_pdf_ui(us[i,], tus[i,], ts_u[i])
    }
    return( lq )
}

################### gradients ###########################################

## gradients for variational alpha prms
grad_tm <- function(alpha, tm, tsig){
    ## calculate the gradient of \log q wrt tilde{m}
    (alpha - tm)/(tsig^2)
}

grad_log_tsig <- function(alpha, tm, tsig){
    ## calculate the gradient of \log q wrt \log tilde{sig} (have >0 constraint)
    - 1 + ( (alpha - tm)^2 / (tsig^2 ) )
}

## gradients for variational coord prms
grad_trstr <- function(ui, tui, tsi){
    ## gradient \log q wrt \tilde{r}*i

    plri = us_to_plr(tui) # convert tui to polar coordinates
    
    ## first calculate dlq_dtphi 
    y_tui = 1 + 2*sum( (tui-ui)^2 ) / ( (1 - sum(ui^2) )*(1 - sum(tui^2) ) )
    dtu_dr = c( cos(plri[2]), sin(plri[2]) )
    ddp_dr = ( y_tui^2 - 1 )^(-1/2) * (4 / (1 - sum(ui^2)) ) * ( sum( dtu_dr*(tui-ui) )/( 1 - sum(tui^2) )  + sum( (tui-ui)^2 )*sum( tui*dtu_dr )/(1 - sum( tui^2) )^2 )
    
    dlq_dtr = - (hypdistpoinc( ui, tui ) / tsi^2) * ddp_dr
    
    ## multiply by dtr_dtrstr
    rstri = phi_to_phistr( plri[1] )
    dlq_dtr_str = dlq_dtr * exp( - rstri ) / ( 1 + exp( - rstri ) )^2

    return( dlq_dtr_str )
}

grad_tphstr <- function(ui, tui, tsi){
    ## gradient \log q wrt \tilde{upvarphi}*i

    plri = us_to_plr(tui) # convert tui to polar coordinates
    
    ## first calculate dlq_dtphi 
    y_tui = 1 + 2*sum( (tui-ui)^2 ) / ( (1 - sum(ui^2) )*(1 - sum(tui^2) ) )
    dtu_dphi = plri[1] * c( -sin(plri[2]), cos(plri[2]) )
    ddp_dphi = ( y_tui^2 - 1 )^(-1/2) * (4 / (1 - sum(ui^2)) ) * ( sum( dtu_dphi*(tui-ui) )/( 1 - sum(tui^2) )  + sum( (tui-ui)^2 )*sum( tui*dtu_dphi )/(1 - sum( tui^2) )^2 )
    
    dlq_tphi = - (hypdistpoinc( ui, tui ) / tsi^2) * ddp_dphi
    
    ## multiply by dtphi_dtphstri
    phstri = phi_to_phistr( plri[2] )
    dlq_dtphstri = dlq_tphi ## * 2 * pi * exp( - phstri ) / ( 1 + exp( - phstri ) )^2

    return( dlq_dtphstri )
}

grad_tphstr_constr <- function(ui, tui, tsi){
    ## gradient \log q wrt \tilde{upvarphi}*i

    plri = us_to_plr(tui) # convert tui to polar coordinates
    
    ## first calculate dlq_dtphi 
    y_tui = 1 + 2*sum( (tui-ui)^2 ) / ( (1 - sum(ui^2) )*(1 - sum(tui^2) ) )
    dtu_dphi = plri[1] * c( -sin(plri[2]), cos(plri[2]) )
    ddp_dphi = ( y_tui^2 - 1 )^(-1/2) * (4 / (1 - sum(ui^2)) ) * ( sum( dtu_dphi*(tui-ui) )/( 1 - sum(tui^2) )  + sum( (tui-ui)^2 )*sum( tui*dtu_dphi )/(1 - sum( tui^2) )^2 )
    
    dlq_tphi = - (hypdistpoinc( ui, tui ) / tsi^2) * ddp_dphi
    
    ## multiply by dtphi_dtphstri
    phstri = phi_to_phistr( plri[2] )
    dlq_dtphstri = dlq_tphi * 2 * pi * exp( - phstri ) / ( 1 + exp( - phstri ) )^2

    return( dlq_dtphstri )
}


grad_log_ts <- function(ui, tui, tsi){
    -1 - tsi^2 + (hypdistpoinc(ui, tui)^2 / tsi^2) - (2*tsi*exp(-tsi^2 /2))/(sqrt(2*pi)*erf(tsi/sqrt(2)))
}

######################## update steps ################################

update_tmu <- function(prm, G, f, h){ 
    ## f and h are 1-dimensional
    ## G = G + mean(f)^2 # scaling with adagrad
    G = .9*G + .1*mean(f)^2 # scaling with rmsprop
    astr = cov(f, h) / var(h)
    eps = (1/ sqrt(G) )*(mean(f - astr*h))
    prm_new = prm + eps
    return( list(prm_new=prm_new, G=G) )
}

update_log_tsig <- function(prm, G, f, h){ 
    G = .9*G + .1*mean(f)^2 # scaling with adagrad
    astr = cov(f, h) / var(h)
    eps = (1/ sqrt(G) )*(mean(f- astr*h))
    prm_new = exp( log( prm ) + eps ) # update on log scale (no >0 constraint)
    return( list(prm_new=prm_new, G=G) )
}

update_tr <- function(prm, G, f, h){ 
    G = .9*G + .1*mean(f)^2 # scaling with adagrad
    astr = cov(f, h) / var(h)
    eps = (1/ sqrt(G) )*(mean(f- astr*h))
    prm_new = rstr_to_r( r_to_rstr( prm ) + eps )
    return( list(prm_new=prm_new, G=G) )
}

update_tphi <- function(prm, G, f, h){ 
    G = .9*G + .1*mean(f)^2 # scaling with adagrad
    astr = cov(f, h) / var(h)
    eps = (1/ sqrt(G) )*(mean(f- astr*h))
    ##prm_new = phistr_to_phi( phi_to_phistr( prm ) + eps )
    prm_new =  prm + eps
    return( list(prm_new=prm_new, G=G) )
}

update_tphi_constr <- function(prm, G, f, h){ 
    G = .9*G + .1*mean(f)^2 # scaling with adagrad
    astr = cov(f, h) / var(h)
    eps = (1/ sqrt(G) )*(mean(f- astr*h))
    prm_new = r_to_pi( pi_to_r( prm ) + eps )
    return( list(prm_new=prm_new, G=G) )
}

############################ initialisation ##############################

apply_isometry <- function(us, pmn, bk_id=c(1,2,3)){
    ## assume constrained coordinates are indexed by (1,2) for now
    ## come back and improve later

    u_bk1 = us[bk_id[1],]
    u_bk2 = us[bk_id[2],]

    ## calculate parameters of transformation
    dist_id = hypdistpoinc(u_bk1, u_bk2)
    a = sqrt( (cosh(dist_id)-1)/(1 + cosh(dist_id)) )

    t1 = u_bk1[1] + 1i*u_bk1[2]
    t2 = a * ( (u_bk1[1] - 1i*u_bk1[2])*(u_bk2[1] + 1i*u_bk2[2]) - 1 ) / ( u_bk2[1] - u_bk1[1] + 1i*(u_bk2[2] - u_bk1[2]) )

    ## apply isometry
    us = t2 * ( us[,1] + 1i*us[,2] - t1 ) / ( ( us[,1] + 1i*us[,2] )*( u_bk1[1] - 1i*u_bk1[2] ) - 1 )

    ## apply isometry to target coordinate mean also
    pmn = t2 * ( pmn[1] + 1i*pmn[2] - t1 ) / ( ( pmn[1] + 1i*pmn[2] )*( u_bk1[1] - 1i*u_bk1[2] ) - 1 )

    ## check 3rd bookstein coordinate positive, reflect if not
    if (Im(us[bk_id[3]]) < 0){
        return( list( us_t = cbind( Re(us), -Im(us) ), pmn_t = c(Re(pmn), -Im(pmn)) ) )
    } else {
        return( list( us_t = cbind( Re(us), Im(us) ), pmn_t = c(Re(pmn), Im(pmn)) ) ) 
    }
}


#################### plotting function #############################

plot_output <- function(out, plt_rng, cols, col_nodes, bk_id=c(1,2,3), labels=NA){
    N = dim(out$tus)[1]
    nd_id = (1:N)[-bk_id[1]]
    maxrep = max(plt_rng) # maximum index - 'final value'
    par(mfrow=c(1,1))
    layout(matrix( c(1,2,3,7,4,5,6,7), nrow=2, byrow=TRUE), widths=c(1,1,1,2) )
    plot(out$tm[plt_rng], type='l', xlab="iteration", ylab=expression(tilde(m)))
    plot(out$tsig[plt_rng], type='l', xlab="iteration", ylab=expression(tilde(sigma)))
    plot(out$ts_u[nd_id[1],plt_rng], type='l', col=cols[1], ylim=c(0, max(out$ts_u)), xlab="iteration", ylab=expression(tilde(s)) )
    for (i in nd_id[-1]){points(out$ts_u[i,plt_rng], type='l', col=cols[i])}
    plot(out$tus_plr[nd_id[1],1,plt_rng], col=col_nodes[1], type='l', ylim=c(0, 1), xlab="iteration", ylab=expression(tilde(r)))
    for(i in nd_id[-1]){points(out$tus_plr[i,1,plt_rng], col=col_nodes[i], type='l')}
    plot(out$tus_plr[nd_id[1],2,plt_rng], col=col_nodes[1], type='l', ylim=c(0, max(out$tus_plr[,2,plt_rng]) ), xlab="iteration", ylab=expression(tilde(phi)))
    for(i in nd_id[-1]){points(out$tus_plr[i,2,plt_rng], col=col_nodes[i], type='l')}
    plot(out$ELBO[plt_rng], type='l', xlab="iteration", ylab="ELBO")

    ## plot latent trajectories, grey to unique colour for each i
    plot( NA, xlim=c(-1,1), ylim=c(-1,1), xlab=expression(ui1), ylab=expression(ui2))
    for (i in 1:N){
        col_nd = colorRampPalette( color = c('grey', col_nodes[i]))
        col_nd = col_nd(maxrep+1)
        points( out$tus[i,1,plt_rng], out$tus[i,2,plt_rng], col=col_nd, pch=19)
    }
    if (is.na(labels)){
        text( out$tus[,1,maxrep], out$tus[,2,maxrep], 1:N, col='black', pos=3)
    }else{
        text( out$tus[,1,maxrep], out$tus[,2,maxrep], labels, col='black', pos=3)
    }

}

####################### black box #####################################

hyperbolic_BBVI <- function(obs, prms_p, S, maxItr, bk_id=c(1,2,3), tol=1e-5){
                                        # obs is a NxN symmetric matrix of observations
                                        # prms_p is list of parameters in posterior
                                        # S in the number of MC samples
                                            # maxItr is maximum number of iterations
                                        # bk_id are constrained coordinates - set a procedure for this?

    N = dim(obs)[1] # number of nodes, obs is NxN adjacency matrix
    
    ## storage for gradient estimation
    f_tm_ts = matrix(0, nrow=S, ncol=2) 
    f_tus = array(0, c(S,N,3)) # have 3 prms for each node (2 angles and k)
    h_tm_ts = matrix(0, nrow=S, ncol=2)
    h_tus = array(0, c(S,N,3))
    G_tm_ts = rep(0,2) # scaling for tm and tsig
    G_tus = matrix(0, nrow=N, ncol=3) # cols: tomega, tphi, tk
    
    ## set up storage for parameters
    hst = init_strg(N, maxItr) 
    
    ## initialise (include isometry transformation for coordinates)
    uinit = init_us_poincare(obs) #prms_p$us
    isom = apply_isometry(uinit, prms_p$mu, bk_id)
    prms_p$mu = isom$pmn_t # 'recenter' target mean
    hst$tus[,,1] = round(isom$us_t,3)
    hst$tus_plr[,,1] = us_to_plr( hst$tus[,,1] )
    hst$ts_u[,1] = rep(0.01, N) # dispersion for each coordinate
    hst$tm[1] = init_alpha_poincare( obs, hst$tus[,,1] )
    hst$tsig[1] = .01
    
    #print(hst$tus[,,1])
    #print(hst$tus_plr[,,1])
    #print(hst$ts_u[,1])
    #print(hst$tm[1])
    #print(hst$tsig[1])
    
    lik_mat<-matrix(0,maxItr,1)
    
    for( irep in 1:maxItr){

        print(irep)
      
        #print(hst$tus[,,irep])
        #print(hst$tm[irep])
        
        lik_mat[irep]<-log_py(Y, hst$tus[,,irep], hst$tm[irep])
      
        #print(lik_mat[irep])
        hst$lik_mat[irep]<-lik_mat[irep]
        ## sample proposals from q
        
       # print(hst$tm[irep])
      #  print(hst$tsig[irep])
      #  print(hst$tus[,,irep])
      #  print(hst$ts_u[,irep])
      #  print(bk_id)
        
        smp = samp_from_q(S, hst$tm[irep], hst$tsig[irep], hst$tus[,,irep], hst$ts_u[,irep], bk_id )
        hst$ELBO[irep] = 0 # set to 0, sum each s and divide by S
       
        
        for (s in 1:S){
            
                        ## evaluate log target
            lp_y = log_py(obs, smp$us_s[,,s], smp$alph_s[s])
            lp_alp = log_post_alpha(smp$alph_s[s], prms_p$m, prms_p$sig)
            lp_us = log_post_us(smp$us_s[,,s], prms_p$mu, prms_p$su)
            ## evaluate log posterior
            lq_alp = log_q_alpha(smp$alph_s[s], hst$tm[irep], hst$tsig[irep])            
            lq_us = log_q_us(smp$us_s[,,s], hst$tus[,,irep], hst$ts_u[,irep], bk_id)

            ## update ELBO
            hst$ELBO[irep] = hst$ELBO[irep] + lp_y + lp_alp + lp_us - lq_alp - lq_us
            
            ## update tm (tilde{m})
            h_tm_ts[s,1] = grad_tm(smp$alph_s[s], hst$tm[irep], hst$tsig[irep])
            f_tm_ts[s,1] = h_tm_ts[s,1] * ( lp_y + lp_alp - lq_alp )
            
            ## update for log tsig (tilde{s})
            h_tm_ts[s,2] = grad_log_tsig(smp$alph_s[s], hst$tm[irep], hst$tsig[irep])
            f_tm_ts[s,2] = h_tm_ts[s,2] * ( lp_y + lp_alp - lq_alp )

            for (i in (1:N)){
               
                ## update for tr_str
                h_tus[s,i,1] = grad_trstr(smp$us_s[i,,s], hst$tus[i,,irep], hst$ts_u[i,irep])
                f_tus[s,i,1] = h_tus[s,i,1] * ( lp_y + lp_us - lq_us )
                
                ## update for tphi_str
                if (i == bk_id[3]){
                    h_tus[s,i,2] = grad_tphstr_constr(smp$us_s[i,,s],hst$tus[i,,irep], hst$ts_u[i,irep])
                    f_tus[s,i,2] = h_tus[s,i,2] * ( lp_y + lp_us - lq_us )
                } else {
                    h_tus[s,i,2] = grad_tphstr(smp$us_s[i,,s], hst$tus[i,,irep], hst$ts_u[i,irep])
                    f_tus[s,i,2] = h_tus[s,i,2] * ( lp_y + lp_us - lq_us )
                }
                
                ## update for log ts_u
                h_tus[s,i,3] = grad_log_ts(smp$us_s[i,,s], hst$tus[i,,irep], hst$ts_u[i,irep])
#print(lp_y)
#print(lp_us)
#print(lq_us)
                f_tus[s,i,3] = h_tus[s,i,3] * ( lp_y + lp_us - lq_us )
            }
        }
        ## ELBO as mean
        hst$ELBO[irep] = hst$ELBO[irep] / S
        
        ## UPDATE STEPS
        tmn_upd = update_tmu(hst$tm[irep], G_tm_ts[1], f_tm_ts[,1], h_tm_ts[,1])
        hst$tm[irep+1] = tmn_upd$prm_new
        G_tm_ts[1] = tmn_upd$G 

        tsig_upd = update_log_tsig(hst$tsig[irep], G_tm_ts[2], f_tm_ts[,2], h_tm_ts[,2])
        hst$tsig[irep+1] = tsig_upd$prm_new
        G_tm_ts[2] = tsig_upd$G

        ## update coordinates
        #nd_upd = (1:N)[-bk_id[1]]
        nd_upd = (1:N)#[-bk_id[1]]
        for (i in nd_upd){
                        
            ## update tri_st and store
            tr_upd = update_tr(hst$tus_plr[i,1,irep], G_tus[i,1], f_tus[,i,1], h_tus[,i,1])
            hst$tus_plr[i,1,irep+1] = tr_upd$prm_new
            G_tus[i,1] = tr_upd$G 
            ## update tphi_str and store
            if (i == bk_id[2]){ # angle not updated for 2nd coordinate
                hst$tus_plr[i,2,irep+1] = hst$tus_plr[i,2,irep]
            } else if (i == bk_id[3]){ # angle constrained for 3rd coordinate
                tph_upd = update_tphi_constr(hst$tus_plr[i,2,irep],G_tus[i,2], f_tus[,i,2], h_tus[,i,2])
                hst$tus_plr[i,2,irep+1] = tph_upd$prm_new
                G_tus[i,2] = tph_upd$G
            } else {
                tph_upd = update_tphi(hst$tus_plr[i,2,irep], G_tus[i,2], f_tus[,i,2], h_tus[,i,2])
                hst$tus_plr[i,2,irep+1] = tph_upd$prm_new
                G_tus[i,2] = tph_upd$G
            }
            ## add in update for us
            hst$tus[i,,irep+1] = plr_to_us(hst$tus_plr[i,,irep+1])
            ## update ts_u and store 
            

            #print( f_tus[,i,3])
            ts_upd = update_log_tsig(hst$ts_u[i,irep], G_tus[i,3], f_tus[,i,3], h_tus[,i,3])
            
            hst$ts_u[i,irep+1] = ts_upd$prm_new
            G_tus[i,3] = ts_upd$G
        }
        
        ## update in hst$
        #hst$tus[-nd_upd,,irep+1] =  hst$tus[-nd_upd,,irep]
        #hst$ts_u[-nd_upd,irep+1] = hst$ts_u[-nd_upd,irep]
        #hst$tus_plr[-nd_upd,,irep+1] = hst$tus_plr[-nd_upd,,irep]

        hst$tus[,,irep+1] =  hst$tus[,,irep]
        hst$ts_u[,irep+1] = hst$ts_u[,irep]
        hst$tus_plr[,,irep+1] = hst$tus_plr[,,irep]
        
        
        ## print to make sure not stuck
        if (irep %% 10 == 0 ){ print(paste("finished updates at rep ", irep) ) }

        ## condition to check whether converged (if so, break)
        ## if ( irep > 2 ){
        ##     if ( hst$ELBO[irep] - hst$ELBO[irep-1] < tol ){
        ##         return( hst )
        ##         break
        ##     }
        ## }
    }

    ## take mod of angles
    hst$tus_plr[,2,] = hst$tus_plr[,2,] %% (2*pi)

    return( hst) # hst stores estimates of parameters
}
