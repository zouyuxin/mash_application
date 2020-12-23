library(mashr)
library(plyr)
estimate_null_correlation = function(data, Ulist, init, max_iter = 30, tol=1,
                                     est_cor = TRUE, track_fit = FALSE, prior = c('nullbiased', 'uniform'), details = FALSE, ...){
  if(class(data) != 'mash'){
    stop('data is not a "mash" object')
  }
  if(!is.null(data$L)){
    stop('We cannot estimate the null correlation for the mash contrast data.')
  }

  prior <- match.arg(prior)
  tracking = list()

  if(missing(init)){
    init = tryCatch(estimate_null_correlation_simple(data, est_cor = est_cor), error = function(e) FALSE)
    if(is.logical(init)){
      warning('Use Identity matrix as the initialize null correlation.')
      init = diag(mashr:::n_conditions(data))
    }
  }

  J = mashr:::n_effects(data)
  m.model = fit_mash_V(data, Ulist, V = init, prior=prior,...)
  pi_s = get_estimated_pi(m.model, dimension = 'all')
  prior.v <- mashr:::set_prior(length(pi_s), prior)

  # compute loglikelihood
  log_liks <- numeric(max_iter+1)
  log_liks[1] <- get_loglik(m.model) #+penalty(prior.v, pi_s)
  V = init

  result = list(V = V, mash.model = m.model)

  niter = 0
  while (niter < max_iter){
    niter = niter + 1
    if(track_fit){
      tracking[[niter]] = result
    }
    # max_V
    V = E_V(data, m.model)/J
    if(est_cor){
      V = cov2cor(V)
    }
    m.model = fit_mash_V(data, Ulist, V, prior=prior, ...)
    pi_s = get_estimated_pi(m.model, dimension = 'all')

    log_liks[niter+1] <- get_loglik(m.model)  #+penalty(prior.v, pi_s)
    delta.ll <- log_liks[niter+1] - log_liks[niter]

    result = list(V = V, mash.model = m.model)
    if (abs(delta.ll) <= tol){
      niter = niter + 1
      break
    }
  }

  log_liks = log_liks[1:niter] #remove tailing NAs
  result$loglik = log_liks
  result$niter = niter
  if(track_fit){
    result$trace = tracking
  }

  if(details){
    return(result)
  }else{
    return(result$V)
  }
}

#' @importFrom plyr aaply laply
E_V = function(data, m.model){
  J = mashr:::n_effects(data)
  Z = data$Bhat/data$Shat
  Shat = data$Shat * data$Shat_alpha
  post.m.shat = m.model$result$PosteriorMean / Shat
  post.sec.shat = laply(1:J, function(i) (t(m.model$result$PosteriorCov[,,i]/Shat[i,])/Shat[i,]) +
                          tcrossprod(post.m.shat[i,])) # JxRxR array
  temp1 = crossprod(Z)
  temp2 = crossprod(post.m.shat, Z) + crossprod(Z, post.m.shat)
  temp3 = unname(aaply(post.sec.shat, c(2,3), sum))

  V = (temp1 - temp2 + temp3)
  # avoid numerical unsymmetry
  V = (V+t(V))/2
}

fit_mash_V <- function(data, Ulist, V, prior=c('nullbiased', 'uniform'), ...){
  data.V = mash_update_data(data, V=V)
  m.model = mash(data.V, Ulist, prior=prior, verbose = FALSE, outputlevel = 3, ...)
  return(m.model)
}
