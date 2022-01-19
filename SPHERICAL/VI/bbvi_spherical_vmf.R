##rm(list= ls()) # removes all from workspace

source("sample_network_functions.R") # function to sample a network
library(igraph)
library(smacof) # for MDS in sphere

## the model is p(y|u,alpha) p(alpha|m,sig) p(u|k,mu) 
## note: coordinates assumed to follow a von-mises fisher distribution (normalising constant has exact form)
## the variational family is q(alpha|tm,tsig) q(u|tk, tmu)

##################### init storage #############################

init_strg <- function(N, maxItr){    
    hst = list(
        ## variational parameters
        tm = rep(NA, maxItr+1), # alpha mean
        tsig = rep(NA, maxItr+1), # alpha variance
        tk = matrix(NA, N, maxItr+1), # spread for coordinates
        tngl = array(NA, c(N,2,maxItr+1)), # polar coordinate means
        tu = array(NA, c(N,3,maxItr+1)), # euclidean coordinate means
        ELBO = rep(NA, maxItr+1) # store ELBO at each iteration
    )
    return( hst )
}

##################### extra useful functions #############################

samp_us_from_vmf <- function(S, tus, tk, bk_id){
    ## sample from a von mises-fisher distn where each coordinate has a unique mean
    ## dimension of sample will be NxdxS  (S is number of samples)

    dim_tus = dim(tus)
    us_smp = array(0, c(dim_tus, S))
    ids_nonbk = (1:dim_tus[1])[-bk_id[1]]
    for (i in ids_nonbk){
        us_smp[i,,] = t(rmovMF( S, tk[i]*tus[i,], alpha = 1))
        ##us_smp[i,,] = t(rvmf( S, mu=tus[i,], k=tk[i]))
    }
    ## keep bk coordinates fixed
    us_smp[bk_id[1],,] = tus[bk_id[1],]
    us_smp[bk_id[2],2,] = tus[bk_id[2],2]  # keep y coord fixed
    
    return( us_smp )
}

samp_from_q <- function(S, tmcur, tsigcur, tuscur, tkcur, bk_id ){
    alph_s = rnorm( S, mean=tmcur, sd=tsigcur ) # alphas
    us_s = samp_us_from_vmf(S, tuscur, tkcur, bk_id )  # coordinates
    return( list(alph_s=alph_s, us_s=us_s) )
}

## functions to convert between cartesian and polar representation of spherical coordinates
us_to_angls <- function(us){
    ## for reference see https://mathworld.wolfram.com/SphericalCoordinates.html
    ## in my notation \theta = phi and \phi = omega
    if (is.vector(us)==TRUE){ us = matrix(us, nrow=1) }
    
    ## calculate angle between z axis
    omega = acos( us[,3] ) 

    ## calculate angle from x axis
    us1pos = us[,1] >= 0
    us2pos = us[,2] >= 0
    # need to add on multiple of pi so that phi \in [0, 2pi) (multiple depends on xy quadrant)
    phi = atan( us[,2]/us[,1] ) + pi*( (!us1pos & us2pos) | (!us1pos & !us2pos) ) + (2*pi)*( us1pos & !us2pos )
    phi[ is.nan(phi) ] = 0 # by convention, take 0 angle (only for c(0,0,1))
    
    return( cbind(omega, phi) )
}

angls_to_us <- function(angles){
    ## ith row of angles is (omega_i, phi_i)
    if (is.vector(angles)==TRUE){ angles = matrix( angles, nrow=1 )}
    us = cbind( sin(angles[,1])*cos(angles[,2]), sin(angles[,1])*sin(angles[,2]), cos(angles[,1]) )
    return(us)
}

## functions to convert between parameters and unconstrained updates (needed for polar angles only)
om_to_omstr <- function( x ){ log( x ) - log( pi - x ) }
omstr_to_om <- function( x ){ pi / (1 + exp( -x )) }
ph_to_phstr <- function( x ){ log( x ) - log( 2*pi - x ) }
phstr_to_ph <- function( x ){ 2*pi / (1 + exp( -x )) }

######################## log densities ##################################

## note: should be easy to vectorise log_py functions
log_py <- function(obs, us, alpha){ 
    ll = 0
    N = dim(us)[1]
    for (i in 1:(N-1)){
        for (j in (i+1):N){
            eta = alpha - acos( sum(us[i,]*us[j,]) )
            ll = ll + obs[i,j]*eta - log(1 + exp(eta))
        }
    }
    
    #print(ll)
    #print(alpha)
    return( ll )
}

log_py_ui <- function(obs, us, alpha, ind){ 
    ## return only the terms in log p(y | u, alpha ) which contain node ind
    ## note: us is s^th sample for node ind, fixed to truth everywhere else    
    ll = 0
    N = dim(us)[1]
    for (j in (1:N)[-ind]){ # loop over all node pairs
            eta = alpha - acos( sum(us[ind,]*us[j,]) )
            ll = ll + obs[ind,j]*eta - log(1 + exp(eta))
    }
    return( ll )
}

log_vMF_pdf <- function(us, tus, tk){
    sum( tk*diag( us %*% t(tus) ) + log(tk) - log(2*pi) - log( exp(tk) - exp(-tk)) )
    ## N = dim(us)[1]
    ## lpdf = 0
    ## for ( i in 1:N){ lpdf = lpdf + log_vMF_pdf_ui(us[i,], mu[i,], k[i]) }
    return( lpdf )
}

log_vMF_pdf_ui <- function(us, mu, k){
    dmovMF( matrix(us, nrow=1), matrix( k*mu, nrow=1), alpha = 1, log = TRUE ) 
    #k*( sum(c(us)*c(mu)) ) + log( k ) - log(2*pi) - log( exp( k ) - exp(- k) )
}

## _post_ functions involve the MODEL PARAMETERS
log_post_alpha <- function(alpha, m, sig){
    ## subset of log posterior which contains alpha:    
    ## log p(alpha | m, sig )
    lp = dnorm(alpha, mean=m, sd=sig, log=TRUE )
    return( lp )
}

log_post_us <- function(us, mu, k){
    ## calculate sum log p(u_i | mu, k ) (same for all i)
    N = dim(us)[1]
    lp = 0 #log_py(obs, us, alpha)
    for ( i in (1:N) ){ lp = lp + log_vMF_pdf_ui(us[i,], mu, k) }
    return( lp )
}

## _q_ functions involve the VARIATIONAL PARAMETERS
log_q_alpha <- function(alpha, tm, tsig){
    ## calculate log q(alpha | tm, tsig )
    lq = dnorm(alpha, mean=tm, sd=tsig, log=TRUE )
    return( lq )
}

log_q_ui <- function(us, tus, tk){
    ## calculate log q( ui | tui, tki )
    ## lq = dmovMF(matrix(us[ind,], nrow=1), theta=matrix(tk[ind]*tus[ind,], nrow=1), alpha=1, log=TRUE)
    lq = log_vMF_pdf_ui(us, tus, tk)
    return( lq )
}

log_q_us <- function(us, tus, tk, bk_id){
    ## calculate log q( ui | tui, tki )
    N = dim(us)[1]
    lq = 0 
    for ( i in (1:N)[-bk_id[1]] ){
        lq = lq + log_vMF_pdf_ui(us[i,], tus[i,], tk[i])
    }
    
    return( lq )
}

################### gradients ###########################################

grad_tm <- function(alpha, tm, tsig){
    ## calculate the gradient of \log q wrt tilde{m}
    (alpha - tm)/(tsig^2)
}

grad_log_tsig <- function(alpha, tm, tsig){
    ## calculate the gradient of \log q wrt \log tilde{sig} (have >0 constraint)
    - 1 + ( (alpha - tm)^2 / (tsig^2 ) )
}

grad_log_tki <- function(ui, tui, tki){
    ## gradient of log q (ui | tui, tki ) wrt log tilde{ki} (have >0 constraint)
    ##1 + tki * sum(c(tui)*c(ui)) - tki* ( exp(tki) + exp(-tki) )/( exp(tki) - exp(-tki) )

    ## this form of the gradient more stable (eg for large tki)
    1 + tki * (sum(c(tui)*c(ui)) - 1) - 2*tki*exp(-2*tki) /( 1 - exp(-2*tki) )
}

grad_tomstr <- function(ui, tki, tangli){
    ## gradient of log q( ui | tui, tki ) with respect to tomega*
    ## tangli is 2d vector = (tomega_i, tphi_i)
    
    dtu_dto = c( cos(tangli[1])*cos(tangli[2]), cos(tangli[1])*sin(tangli[2]), - sin(tangli[1]) )
    ## not using constraints on angles...
    ## omegstr = om_to_omstr( tangli[1] )
    ## dto_dtos = (pi * exp(- omegstr)) / (1 + exp(- omegstr))^2
    grad = tki * ( dtu_dto %*% ui ) #* dto_dtos # using chain rule 
    
    return( grad )
}

grad_tphstr <- function(ui, tki, tangli){
    ## gradient of log q( ui | tui, tki ) with respect to tphi*
    ## tangli is 2d vector = (tomegai, tphii)
    
    dtu_dtph = c( -sin(tangli[1])*sin(tangli[2]), sin(tangli[1])*cos(tangli[2]), 0 )
    ## not using constraints on angles...
    ## phistr = ph_to_phstr( tangli[2] )
    ## dtph_dtphs = ( 2 * pi * exp(- phistr) ) / (1 + exp(- phistr))^2
    grad = tki * ( dtu_dtph %*% ui ) #* dtph_dtphs # using chain rule here

    return( grad )
}


######################## update steps ################################

update_tmn <- function(tmcur, G, f, h){
    ## f and h are 1-dimensional
    ##G = G + mean(f)^2 # scaling with adagrad
    G = .9*G + .1*mean(f)^2 # scaling with rmsprop
    astr = cov(f, h) / var(h) 
    eps = (.5/ sqrt(G) )*(mean(f - astr*h))
    tm_new = tmcur + eps
    return( list(tm_new=tm_new, G=G) )
}

update_tsig <- function(tsigcur, G, f, h){
    ##G = G + mean(f)^2 # scaling with adagrad
    G = .9*G + .1*mean(f)^2 # scaling with rmsprop
    astr = cov(f, h) / var(h) 
    eps = (.5/ sqrt(G) )*(mean(f- astr*h))
    tsig_new = exp( log(tsigcur) + eps ) # update on log scale (no >0 constraint)
    return( list(tsig_new=tsig_new, G=G) )
}

update_tom <- function(tomcur, G, f, h){
    ##G = G + mean(f)^2
    G = .9*G + .1*mean(f)^2 # scaling with rmsprop
    astr = cov(f, h) / var(h)
    if (is.nan(astr)){ astr = 1 }
    eps = (.75/ sqrt(G) )*(mean(f- astr*h))
    ##tom_new = omstr_to_om( om_to_omstr(tomcur) + eps )
    tom_new = tomcur + eps 
    return( list(tom_new=tom_new, G=G) )
}

update_tph <- function(tphcur, G, f, h){
    ##G = G + mean(f)^2
    G = .9*G + .1*mean(f)^2 # scaling with rmsprop
    astr = cov(f, h) / var(h)
    eps = (.75/ sqrt(G) )*(mean(f - astr*h))
    ##tph_new = phstr_to_ph( ph_to_phstr(tphcur) + eps )
    tph_new = tphcur + eps
    return( list(tph_new=tph_new, G=G) )
}

update_tk <- function(tkcur, G, f, h){
    ##G = G + mean(f)^2 # scaling with adagrad
    G = .9*G + .1*mean(f)^2 # scaling with rmsprop
    astr = cov(f, h) / var(h)
    if (is.nan(astr) & (cov(f,h)==var(h))){ astr = 1 }
    eps = (.5/ sqrt(G) )*(mean(f- astr*h))
    tk_new = exp( log(tkcur) + eps ) # update on log scale (no >0 constraint)
    return( list(tk_new=tk_new, G=G) )
}

##################### initialisation #############################

rotmat <- function( thtx, thty, thtz ){
    matrix( c( cos(thtz)*cos(thty), sin(thtz)*cos(thty), - sin(thty), cos(thtz)*sin(thty)*sin(thtx) - sin(thtz)*cos(thtx), sin(thtz)*sin(thty)*sin(thtx) + cos(thtz)*cos(thtx), cos(thty)*sin(thtx), cos(thtz)*sin(thty)*cos(thtx) + sin(thtz)*sin(thtx), sin(thtz)*sin(thty)*cos(thtx) - cos(thtz)*sin(thtx), cos(thty)*cos(thtx)), ncol=3, byrow=FALSE )
}

apply_isometry <- function(us, pmn, bk_id = c(1,2)){
    
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

    ## apply transform to mean of target
    pmn_t = c( R %*% t(pmn) )
    
    return( list(us_t = us_t, pmn_t=pmn_t) )
}

init_us_sphere <- function(obs){
    ## initialise coordinates for spherical coordinates, assumes 3d

    ## get graph distances
    g_obs = graph_from_adjacency_matrix( obs, mode="undirected" )
    ds = distances( g_obs )
    
    ## remove Inf if disconnected
    mx = max( ds[ ds < Inf ] ) # largest finite observed value 
    ds[ ds==Inf ] = mx + 1 
   
    ## get mds coords 
    us = smacofSphere( ds, 3, type="ratio")$conf
    us = us / sqrt(rowSums(us^2))  # to ensure lie exactly on sphere
    
    return( us )    
}

init_alpha_sphere <- function( obs, us_init, n_alph=20 ){
    ## note: assumes 3-dimensional
    
    ## simple grid search, take the most likely
    ll = rep(0, n_alph)
    avec = seq(0, 3, length=n_alph)
    for (i in 1:n_alph){ ll[i] = log_py(obs, us_init, avec[i]) }

    a_init = avec[ which.max(ll) ] # take 'most likely' value

    return( a_init )
}

####################### plot output #####################################

plot_traces <- function(output, plt_rng, cols, col_nodes, bk_id=c(1,2)){
    N = dim(output$tk)[1]
    nd_id = (1:N)[-bk_id[1]]
    ## plot variational parameters
    par(mfrow=c(2,3))
    plot(output$tm[plt_rng], type='l', xlab="iteration", ylab=expression(tilde(m)) )
    plot(output$tsig[plt_rng], type='l', xlab="iteration", ylab=expression(tilde(sigma)) )
    plot(output$tk[nd_id[1],plt_rng], type='l', col=cols[1], ylim=c(0, max(output$tk)), xlab="iteration", ylab=expression(tilde(k)))
    for (i in nd_id[-1]){points(output$tk[i,plt_rng], type='l', col=cols[i])}
    plot(output$tngl[nd_id[1],1,plt_rng], col=col_nodes[1], type='l', ylim=c(0, pi), xlab="iteration", ylab=expression(tilde(phi)))
    for(i in nd_id[-1]){points(output$tngl[i,1,plt_rng], col=col_nodes[i], type='l')}
    plot(output$tngl[nd_id[1],2,plt_rng], col=col_nodes[1], type='l', ylim=c(0, 2*pi), xlab="iteration", ylab=expression(tilde(omega)))
    for(i in nd_id[-1]){points(output$tngl[i,2,plt_rng], col=col_nodes[i], type='l')}
    plot(output$ELBO[plt_rng], type='l', xlab="iteration", ylab="ELBO")
}

plot_spherical_coords <- function(tuTrc, plt_rng, col_nodes, labels=NA){
    N = dim(tuTrc)[1]
    maxrep = max(plt_rng) # largest index
    clear3d()  
    spheres3d(0,0,0,lit=FALSE,color="white")
    spheres3d(0,0,0,radius=1.01,lit=FALSE,color="black",front="lines")
    for (i in 1:N){
        col_nd = colorRampPalette( color = c('grey', col_nodes[i]))
        col_nd = col_nd(maxrep+1)
        spheres3d(tuTrc[i,1,plt_rng], tuTrc[i,2,plt_rng], tuTrc[i,3,plt_rng], col=col_nd, radius=0.03)
    }
    if (is.na(labels)){
        text3d(tuTrc[,1,maxrep], tuTrc[,2,maxrep], tuTrc[,3,maxrep], 1:N, pos=3)
    } else{
        text3d(tuTrc[,1,maxrep], tuTrc[,2,maxrep], tuTrc[,3,maxrep], labels, pos=3)
    }
}

####################### black box EM #####################################

spherical_BBVI <- function(obs, prms_p, S, maxItr, bk_id=c(1,2), tol=1e-5){
                                        # obs is a NxN symmetric matrix of observations
                                        # prms_p is list of parameters in posterior (to be chosen by the user)
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
    
    ## initialise : first start everything at the truth
    ## u_init = prms_p$us 
    u_init = init_us_sphere( obs )
    prms_p$mu = movMF( u_init, 1)$theta # estimate this from init (smacof can center anywhere...)
    prms_p$mu = prms_p$mu / sqrt( sum(prms_p$mu^2) )
    a_init = init_alpha_sphere( obs, u_init, n_alph=20 )
    isom = apply_isometry(u_init, prms_p$mu, bk_id)
    prms_p$mu = isom$pmn_t # update 'center' of coords
    hst$tm[1] = a_init # max ll init
    hst$tsig[1] = .1
    hst$tk[,1] = rep(10,N) # reasonably concentrated
    hst$tu[,,1] = isom$us_t # transformed MDS init
    hst$tngl[,,1] = us_to_angls( hst$tu[,,1] ) 
    
    for( irep in 1:maxItr){

        ## ##################################
        ## update tilde parameters
        ## ##################################
        
        ## sample proposals from q 
        smp = samp_from_q(S, hst$tm[irep], hst$tsig[irep], hst$tu[,,irep], hst$tk[,irep], bk_id)
        hst$ELBO[irep] = 0 # set to 0, sum in loop and then divide by S
        
        for (s in 1:S){
            ## evaluate log target
            lp_y = log_py(obs, smp$us_s[,,s], smp$alph_s[s])
            lp_us = log_post_us(smp$us_s[,,s], prms_p$mu, prms_p$k)
            lp_alp = log_post_alpha(smp$alph_s[s], prms_p$m, prms_p$sig)
            ## evaluate log variational family
            lq_alp = log_q_alpha(smp$alph_s[s], hst$tm[irep], hst$tsig[irep])
            lq_us = log_q_us(smp$us_s[,,s], hst$tu[,,irep], hst$tk[,irep], bk_id)

            ## calculate ELBO
            hst$ELBO[irep] = hst$ELBO[irep] + lp_y + lp_us + lp_alp - lq_alp - lq_us
            
            ## tm grad
            h_tm_ts[s,1] = grad_tm(smp$alph_s[s], hst$tm[irep], hst$tsig[irep]) 
            f_tm_ts[s,1] = h_tm_ts[s,1] * ( lp_y + lp_alp - lq_alp )

            ## log tsig grad (tilde{s})
            h_tm_ts[s,2] = grad_log_tsig(smp$alph_s[s], hst$tm[irep], hst$tsig[irep])
            f_tm_ts[s,2] = h_tm_ts[s,2] * ( lp_y + lp_alp - lq_alp )

            for (i in (1:N)){ 
                ## update for tomeg_str
                h_tus[s,i,1] = grad_tomstr(smp$us_s[i,,s], hst$tk[i,irep], hst$tngl[i,,irep])
                f_tus[s,i,1] = h_tus[s,i,1] * ( lp_y + lp_us - lq_us )

                ## update for tphi_str
                h_tus[s,i,2] = grad_tphstr(smp$us_s[i,,s], hst$tk[i,irep], hst$tngl[i,,irep])
                f_tus[s,i,2] = h_tus[s,i,2] * ( lp_y + lp_us - lq_us )
                
                ## update for log tki (tilde{ki})
                h_tus[s,i,3] = grad_log_tki(smp$us_s[i,,s], hst$tu[i,,irep], hst$tk[i,irep])
                f_tus[s,i,3] = h_tus[s,i,3] * ( lp_y + lp_us - lq_us )
            }
        }
        ## take ELBO as mean
        hst$ELBO[irep] = hst$ELBO[irep]/S
        
        ## 1) update tm and tsig
        tmn_upd = update_tmn(hst$tm[irep], G_tm_ts[1], f_tm_ts[,1], h_tm_ts[,1])
        hst$tm[irep+1] = tmn_upd$tm_new
        G_tm_ts[1] = tmn_upd$G 
        tsig_upd = update_tsig(hst$tsig[irep], G_tm_ts[2], f_tm_ts[,2], h_tm_ts[,2])
        hst$tsig[irep+1] = tsig_upd$tsig_new
        G_tm_ts[2] = tsig_upd$G       

        nd_upd = (1:N)[-bk_id[1]]
        for (i in nd_upd){       
            ## update tomeg_st and store
            tom_upd = update_tom(hst$tngl[i,1,irep], G_tus[i,1], f_tus[,i,1], h_tus[,i,1])
            hst$tngl[i,1,irep+1] = tom_upd$tom_new
            ##hst$tngl[i,1,irep+1] = hst$tngl[i,1,irep]
            G_tus[i,1] = tom_upd$G             
            ## update tphi_str and store
            if (i == bk_id[2]){
                ## tilde phi constant
                hst$tngl[i,2,irep+1] = hst$tngl[i,2,irep] 
            } else {
                tph_upd = update_tph(hst$tngl[i,2,irep], G_tus[i,2], f_tus[,i,2], h_tus[,i,2])
                hst$tngl[i,2,irep+1] = tph_upd$tph_new
                ##hst$tngl[i,2,irep+1] = hst$tngl[i,2,irep]
                G_tus[i,2] = tph_upd$G 
            }
            ## add in update for coordinates (not in polar)
            hst$tu[i,,irep+1] = angls_to_us(hst$tngl[i,,irep+1]) 
            ## update tk and store
            tk_upd = update_tk(hst$tk[i,irep], G_tus[i,3], f_tus[,i,3], h_tus[,i,3])
            hst$tk[i,irep+1] = tk_upd$tk_new
            ##hst$tk[i,irep+1] = hst$tk[i,irep]
            G_tus[i,3] = tk_upd$G
        }

        ## update bk coordinate parameters
        hst$tk[-nd_upd,irep+1] = hst$tk[-nd_upd,irep]
        hst$tu[-nd_upd,,irep+1] = hst$tu[-nd_upd,,irep]
        hst$tngl[-nd_upd,,irep+1] = hst$tngl[-nd_upd,,irep]
        
        ## print to make sure not stuck
        if (irep %% 50 == 0 ){ print(paste("finished updates at rep ", irep) ) }

        ## condition to check whether converged (if so, break)
        ## if ( irep >= 2 ){
        ##     if ( hst$ELBO[irep] - hst$ELBO[irep-1] < tol ){
        ##         return( list(hst=hst, prms_p=prms_p) )
        ##         break
        ##     }
        ## }
    }

    ## apply mod to angles
    hst$tngl[,1,] = hst$tngl[,1,] %% pi
    hst$tngl[,2,] = hst$tngl[,2,] %% (2*pi)
    
    return( list(hst=hst, prms_p=prms_p) ) # hst stores estimates of parameters
}

