simseq <- function(N, prior.mu=0, prior.sd, wt, mucut=0, mucutf=0.05,
                postcut=0.95, postcutf=0.9,
                ignore=20, nsim=1000) {
  prior.mu <- rep(prior.mu, length=2)
  prior.sd <- rep(prior.sd, length=2)
  sd1 <- prior.sd[1]; sd2 <- prior.sd[2]
  v1 <- sd1 ^ 2
  v2 <- sd2 ^ 2
  j <- 1 : N
  cmean <- Mu <- PostN <- Post <- Postf <- postfe <- postmean <- numeric(nsim)
  stopped <- stoppedi <- stoppedf <- stoppedfu <- stopfe <- status <-
    integer(nsim)
  notignored <- - (1 : ignore)

  # Derive function to compute posterior mean
  pmean <- gbayesMixPost(NA, NA, d0=prior.mu[1], d1=prior.mu[2],
                                 v0=v1, v1=v2, mix=wt, what='postmean')
  
  for(i in 1 : nsim) {
    # See http://stats.stackexchange.com/questions/70855
    component <- if(wt == 1) 1 else sample(1 : 2, size=1, prob=c(wt, 1. - wt))
    mu <- prior.mu[component] + rnorm(1) * prior.sd[component]
    # mu <- rnorm(1, mean=prior.mu, sd=prior.sd) if only 1 component
    
    Mu[i] <- mu
    y  <- rnorm(N, mean=mu, sd=1)
    ybar <- cumsum(y) / j    # all N means for N sequential analyses
    pcdf <- gbayesMixPost(ybar, 1. / j,
                          d0=prior.mu[1], d1=prior.mu[2],
                          v0=v1, v1=v2, mix=wt, what='cdf')
    post  <- 1 - pcdf(mucut)
    PostN[i] <- post[N]
    postf <- pcdf(mucutf)
    s <- stopped[i] <-
      if(max(post) < postcut) N else min(which(post >= postcut))
    Post[i]  <- post[s]   # posterior at stopping
    cmean[i] <- ybar[s]   # observed mean at stopping
    # If want to compute posterior median at stopping:
    #    pcdfs <- pcdf(mseq, x=ybar[s], v=1. / s)
    #    postmed[i] <- approx(pcdfs, mseq, xout=0.5, rule=2)$y
    #    if(abs(postmed[i]) == max(mseq)) stop(paste('program error', i))
    postmean[i] <- pmean(x=ybar[s], v=1. / s)
    
    # Compute stopping time if ignore the first "ignore" looks
    stoppedi[i] <- if(max(post[notignored]) < postcut) N
    else
      ignore + min(which(post[notignored] >= postcut))
    
    # Compute stopping time if also allow to stop for futility:
    # posterior probability mu < 0.05 > 0.9
    stoppedf[i] <- if(max(post) < postcut & max(postf) < postcutf) N
    else
      min(which(post >= postcut | postf >= postcutf))
    
    # Compute stopping time for pure futility analysis
    s <- if(max(postf) < postcutf) N else min(which(postf >= postcutf))
    Postf[i] <- postf[s]
    stoppedfu[i] <- s

    ## Another way to do this: find first look that stopped for either
    ## efficacy or futility.  Record status: 0:not stopped early,
    ## 1:stopped early for futility, 2:stopped early for efficacy
    ## Stopping time: stopfe, post prob at stop: postfe

    stp <- post >= postcut | postf >= postcutf
    s <- stopfe[i] <- if(any(stp)) min(which(stp)) else N
    status[i] <- if(any(stp)) ifelse(postf[s] >= postcutf, 1, 2) else 0
    postfe[i] <- if(any(stp)) ifelse(status[i] == 2, post[s],
                                     postf[s]) else post[N]
  }
  list(mu=Mu, post=Post, postn=PostN, postf=Postf,
       stopped=stopped, stoppedi=stoppedi,
       stoppedf=stoppedf, stoppedfu=stoppedfu,
       cmean=cmean, postmean=postmean,
       postfe=postfe, status=status, stopfe=stopfe)
}
