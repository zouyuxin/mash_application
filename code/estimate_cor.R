library(plyr)

#' @param rho the off diagonal element of V, 2 by 2 correlation matrix
#' @param Ulist a list of covariance matrices, U_{k}
get_sigma <- function(rho, Ulist){
  V <- matrix(c(1,rho,rho,1), 2,2)
  lapply(Ulist, function(U) U + V)
}

penalty <- function(prior, pi_s){
  subset <- (prior != 1.0)
  sum((prior-1)[subset]*log(pi_s[subset]))
}

#' @title compute log likelihood
#' @param L log likelihoods,
#' where the (i,k)th entry is the log probability of observation i
#' given it came from component k of g
#' @param p the vector of mixture proportions
#' @param prior the weight for the penalty
compute.log.lik <- function(lL, p, prior){
  p = normalize(pmax(0,p))
  temp = log(exp(lL$loglik_matrix) %*% p)+lL$lfactors
  return(sum(temp) + penalty(prior, p))
  # return(sum(temp))
}

normalize <- function(x){
  x/sum(x)
}
# ---------------------------------------- #
# Optim
# ---------------------------------------- #
#' @title Optimize rho with several initial values
#' @param X data, Z scores
#' @param Ulist a list of covariance matrices (expand)
#' @param init_rho initial value for rho. The user could provide several initial values as a vector.
#' @param prior indicates what penalty to use on the likelihood, if any
#' @return list of result
#' \item{result}{result from the rho which gives the highest log likelihood}
#' \item{status}{whether the result is global max or local max}
#' \item{loglik}{the loglikelihood value}
#' \item{rho}{the estimated rho}
#' \item{time}{the running time for each initial rho}
#'
optimize_pi_rho_times <- function(X, Ulist, init_rho=0, prior=c("nullbiased", "uniform"), tol=1e-5){
  times = length(init_rho)
  result = list()
  loglik = c()
  rho = c()
  time.t = c()
  for(i in 1:times){
    out.time = system.time(result[[i]] <- optimize_pi_rho(X, Ulist,
                                                          init_rho=init_rho[i],
                                                          prior=prior,
                                                          tol=tol))
    time.t = c(time.t, out.time['elapsed'])
    loglik = c(loglik, tail(result[[i]]$loglik, n=1))
    rho = c(rho, result[[i]]$rho)
  }
  if(abs(max(loglik) - min(loglik)) < 1e-4){
    status = 'global'
  }else{
    status = 'local'
  }
  ind = which.max(loglik)
  return(list(result = result[[ind]], status = status, loglik = loglik, time = time.t, rho=rho))
}

#' @title optimize rho
#' @param X data, Z scores
#' @param Ulist a list of covariance matrices
#' @param init_rho an initial value for rho
#' @param tol tolerance for optimizaiton stop
#' @param prior indicates what penalty to use on the likelihood, if any
#' @return list of result
#' \item{pi}{estimated pi}
#' \item{rho}{estimated rho}
#' \item{loglik}{the loglikelihood value at each iteration}
#' \item{niter}{the number of iteration}
#'
optimize_pi_rho <- function(X, Ulist, init_rho=0, tol=1e-5, prior=c("nullbiased", "uniform")){
  prior <- match.arg(prior)
  if(length(Ulist) <= 1){
    stop('Please provide more U! With only one U, the correlation could be estimated directly using mle.')
  }
  prior <- mashr:::set_prior(length(Ulist), prior)

  Sigma <- get_sigma(init_rho, Ulist)
  lL <- t(plyr::laply(Sigma,function(U){mvtnorm::dmvnorm(x=X,sigma=U, log=TRUE)}))
  lfactors    <- apply(lL,1,max)
  matrix_llik <- lL - lfactors
  lL = list(loglik_matrix = matrix_llik,
            lfactors   = lfactors)
  pi_s <- mashr:::optimize_pi(exp(lL$loglik_matrix),prior=prior,optmethod='mixIP')

  log_liks <- c()
  ll       <- compute.log.lik(lL, pi_s, prior)
  log_liks <- c(log_liks, ll)
  delta.ll <- 1
  niter <- 0
  rho_s <- init_rho
  while( delta.ll > tol){
    # max_rho
    rho_s <- optim(rho_s, optimize_rho, lower = -1, upper = 1, X = X, Ulist=Ulist, pi_s = pi_s, prior = prior, method = 'Brent')$par

    Sigma <- get_sigma(rho_s, Ulist)
    lL <- t(plyr::laply(Sigma,function(U){mvtnorm::dmvnorm(x=X,sigma=U, log=TRUE)}))
    lfactors <- apply(lL,1,max)
    matrix_llik <- lL - lfactors
    lL = list(loglik_matrix = matrix_llik,
              lfactors   = lfactors)

    # max pi
    pi_s <- mashr:::optimize_pi(exp(lL$loglik_matrix),prior=prior,optmethod='mixIP')

    # compute loglike
    ll <- compute.log.lik(lL, pi_s, prior)
    log_liks <- c(log_liks, ll)
    # Update delta
    delta.ll <- log_liks[length(log_liks)] - log_liks[length(log_liks)-1]
    niter <- niter + 1
  }
  return(list(pi = pi_s, rho=rho_s, loglik = log_liks, niter = niter))
}

optimize_rho <- function(rho, X, Ulist, pi_s, prior){
  Sigma <- get_sigma(rho, Ulist)
  lL <- t(plyr::laply(Sigma,function(U){mvtnorm::dmvnorm(x=X,sigma=U, log=TRUE)}))
  lfactors <- apply(lL,1,max)
  matrix_llik <- lL - lfactors
  lL = list(loglik_matrix = matrix_llik,
            lfactors   = lfactors)

  return(-compute.log.lik(lL, pi_s, prior))
}

# ---------------------------#
# EM rho optim
# ---------------------------#
mixture.EM.times <- function(X, Ulist, init_rho=0, init_pi=NULL, prior = c('nullbiased', 'uniform'), control = list()){
  times = length(init_rho)
  result = list()
  loglik = c()
  rho = c()
  time.t = c()
  converge.status = c()
  for(i in 1:times){
    out.time = system.time(result[[i]] <- mixture.EM(X, Ulist,
                                                     init_pi=init_pi,
                                                     init_rho=init_rho[i],
                                                     prior=prior,
                                                     control = control))
    time.t = c(time.t, out.time['elapsed'])
    loglik = c(loglik, result[[i]]$loglik)
    rho = c(rho, result[[i]]$rho)
    converge.status = c(converge.status, result[[i]]$converged)
  }
  if(abs(max(loglik) - min(loglik)) < 1e-4){
    status = 'global'
  }else{
    status = 'local'
  }
  ind = which.max(loglik)
  return(list(result = result[[ind]], status = status, loglik = loglik, rho=rho, time = time.t, converge.status = converge.status))
}

mixture.EM <- function(X, Ulist, init_rho=0, init_pi = NULL, prior = c('nullbiased', 'uniform'), control = list()) {
  prior = match.arg(prior)
  prior <- mashr:::set_prior(length(Ulist), prior)
  k = length(Ulist)
  if (is.null(init_pi)){
    init_pi <- rep(1/k,k)
  }
  control = ashr:::set_control_squarem(control,nrow(X))
  res = SQUAREM::squarem(par=c(init_pi, init_rho),fixptfn=fixpoint_EM, objfn=negpenloglik,X=X, Ulist=Ulist, prior=prior, control=control)

  return(list(pihat = normalize(pmax(0,head(res$par, -1))), rho = tail(res$par, 1), loglik=-res$value.objfn, niter = res$iter, converged=res$convergence, control=control))
}

fixpoint_EM = function(par, X, Ulist, prior){
  rho = tail(par,1)
  pi_s = head(par, -1)
  pi_s = normalize(pmax(0,pi_s)) #avoid occasional problems with negative pis due to rounding

  # compute L
  Sigma <- get_sigma(rho, Ulist)
  L <- t(plyr::laply(Sigma,function(U){mvtnorm::dmvnorm(x=X,sigma=U)}))

  # E
  m  = t(pi_s * t(L)) # matrix_lik is n by k; so this is also n by k
  m.rowsum = rowSums(m)
  classprob = m/m.rowsum #an n by k matrix

  # M
  pinew = normalize(colSums(classprob) + prior - 1)

  rhonew = optimize(EMloglikelihood, interval = c(-1,1), maximum = TRUE, X = X, Ulist = Ulist, z = classprob)$maximum

  return(c(pinew,rhonew))
}

EMloglikelihood = function(rho, X, Ulist, z){
  Sigma = get_sigma(rho, Ulist)
  L = t(plyr::laply(Sigma,function(U){mvtnorm::dmvnorm(x=X,sigma=U, log=TRUE)}))
  sum(L * z)
}

negpenloglik = function(par, X, Ulist, prior){
  Sigma <- get_sigma(tail(par,1), Ulist)
  lL <- t(plyr::laply(Sigma,function(U){mvtnorm::dmvnorm(x=X,sigma=U, log=TRUE)}))
  lfactors <- apply(lL,1,max)
  matrix_llik <- lL - lfactors
  lL = list(loglik_matrix = matrix_llik,
            lfactors   = lfactors)
  ll <- compute.log.lik(lL, head(par, -1), prior)
  return(-ll)
}

# ---------------------------------------- #
# EM M rho
# ---------------------------------------- #

mixture.M.rho.times <- function(X, Ulist, init_rho=0, tol=1e-5, prior = c('nullbiased', 'uniform')){
  times = length(init_rho)
  result = list()
  loglik = c()
  rho = c()
  time.t = c()
  converge.status = c()
  for(i in 1:times){
    out.time = system.time(result[[i]] <- mixture.M.rho(X, Ulist,
                                                      init_rho=init_rho[i],
                                                      prior=prior,
                                                      tol=tol))
    time.t = c(time.t, out.time['elapsed'])
    loglik = c(loglik, tail(result[[i]]$loglik, 1))
    rho = c(rho, result[[i]]$rho)
  }
  if(abs(max(loglik) - min(loglik)) < 1e-4){
    status = 'global'
  }else{
    status = 'local'
  }
  ind = which.max(loglik)
  return(list(result = result[[ind]], status = status, loglik = loglik, rho=rho, time = time.t))
}

mixture.M.rho <- function(X, Ulist, init_rho=0, tol=1e-5, prior = c('nullbiased', 'uniform')) {

  prior <- match.arg(prior)

  m.model = fit_mash(X, Ulist, rho = init_rho, prior=prior)
  pi_s = get_estimated_pi(m.model, dimension = 'all')
  prior.v <- mashr:::set_prior(length(pi_s), prior)

  # compute loglikelihood
  loglik <- c()
  loglik <- c(loglik, get_loglik(m.model)+penalty(prior.v, pi_s))
  delta.ll <- 1
  niter <- 0
  rho = init_rho

  while(delta.ll > tol){
    # max_rho
    rho <- E_rho(X, m.model)

    m.model = fit_mash(X, Ulist, rho, prior=prior)

    pi_s = get_estimated_pi(m.model, dimension = 'all')

    loglik <- c(loglik, get_loglik(m.model)+penalty(prior.v, pi_s))
    # Update delta
    delta.ll <- loglik[length(loglik)] - loglik[length(loglik)-1]
    niter <- niter + 1
  }

  return(list(pihat = normalize(pi_s), rho = rho, loglik=loglik))
}

E_rho <- function(X, m.model){
  n = nrow(X)
  post.m = m.model$result$PosteriorMean
  post.sec = plyr::laply(1:n, function(i) m.model$result$PosteriorCov[,,i] + tcrossprod(post.m[i,])) # nx2x2 array

  temp2 = -sum(X[,1]*X[,2]) + sum(X[,1]*post.m[,2]) + sum(X[,2]*post.m[,1]) - sum(post.sec[,1,2])
  temp1 = sum(X[,1]^2 + X[,2]^2) - 2*sum(X[,1]*post.m[,1]) - 2*sum(X[,2]*post.m[,2]) + sum(post.sec[,1,1] + post.sec[,2,2])

  rts = polyroot(c(temp2, temp1-n, temp2, n))

  # check complex number
  is.real = abs(Im(rts))<1e-12
  if(sum(is.real) == 1){
    return(Re(rts[is.real]))
  }else{
    print('3 real roots')
    return(Re(rts))
  }
}

fit_mash <- function(X, Ulist, rho, prior=c('nullbiased', 'uniform')){
  m.data = mashr::mash_set_data(Bhat=X, Shat=1, V = matrix(c(1, rho, rho, 1), 2, 2))
  m.model = mashr::mash(m.data, Ulist, prior=prior, verbose = FALSE, outputlevel = 3)
  return(m.model)
}

# ---------------------------------------- #
# EM V
# ---------------------------------------- #
mixture.MV.times <- function(X, Ulist, init_V = list(diag(ncol(X))), tol=1e-5, prior = c('nullbiased', 'uniform')){
  times = length(init_V)
  result = list()
  loglik = c()
  V = list()
  time.t = c()
  converge.status = c()
  for(i in 1:times){
    out.time = system.time(result[[i]] <- mixture.MV(X, Ulist,
                                                      init_V=init_V[[i]],
                                                      prior=prior,
                                                      tol = tol))
    time.t = c(time.t, out.time['elapsed'])
    loglik = c(loglik, tail(result[[i]]$loglik, 1))
    V = c(V, list(result[[i]]$V))
  }
  if(abs(max(loglik) - min(loglik)) < 1e-4){
    status = 'global'
  }else{
    status = 'local'
  }
  ind = which.max(loglik)
  return(list(result=result[[ind]], status = status, loglik = loglik, V=V, time = time.t))
}

mixture.MV <- function(X, Ulist, init_V=diag(ncol(X)), tol=1e-5, prior = c('nullbiased', 'uniform')){
  prior <- match.arg(prior)

  m.model = fit_mash_V(X, Ulist, V = init_V, prior=prior)
  pi_s = get_estimated_pi(m.model, dimension = 'all')
  prior.v <- mashr:::set_prior(length(pi_s), prior)

  # compute loglikelihood
  log_liks <- c()
  log_liks <- c(log_liks, get_loglik(m.model)+penalty(prior.v, pi_s))
  delta.ll <- 1
  niter <- 0
  V = init_V

  while(delta.ll > tol){
    # max_V
    V = E_V(X, m.model)
    V = cov2cor(V)
    m.model = fit_mash_V(X, Ulist, V, prior=prior)
    pi_s = get_estimated_pi(m.model, dimension = 'all')

    log_liks <- c(log_liks, get_loglik(m.model)+penalty(prior.v, pi_s))
    # Update delta
    delta.ll <- log_liks[length(log_liks)] - log_liks[length(log_liks)-1]
    niter <- niter + 1
  }
  return(list(pi = normalize(pi_s), V=V, loglik = log_liks))
}

E_V = function(X, m.model){
  n = nrow(X)
  post.m = m.model$result$PosteriorMean
  post.sec = plyr::laply(1:n, function(i) m.model$result$PosteriorCov[,,i] + tcrossprod(post.m[i,])) # nx2x2 array

  temp1 = crossprod(X)
  temp2 = crossprod(post.m, X) + crossprod(X, post.m)
  temp3 = unname(aaply(post.sec, c(2,3), sum))

  (temp1 - temp2 + temp3)/n
}

fit_mash_V <- function(X, Ulist, V, prior=c('nullbiased', 'uniform')){
  m.data = mashr::mash_set_data(Bhat=X, Shat=1, V = V)
  m.model = mashr::mash(m.data, Ulist, prior=prior, verbose = FALSE, outputlevel = 3)
  return(m.model)
}

# ---------------------------------------- #
# Estimate Correlation matrix V
# ---------------------------------------- #

#' @title Estimate correlation matrix
#' @description Estimate the pairwise correlation using mle
#' @param data a mash data object
#' @param Ulist a list of covariance matrices
#' @param gridmult scalar indicating factor by which adjacent grid values should differ; close to 1 for fine grid
#' @param grid vector of grid values to use (scaling factors omega in paper)
#' @param normalizeU whether or not to normalize the U covariances to have maximum of 1 on diagonal
#' @param usepointmass whether to include a point mass at 0, corresponding to null in every condition
#' @param init_pi a vector of inital value for pi, mixture proportions (only one initial value)
#' @param init_rho initial value for rho. The user could provide several initial values as a vector.
#' @param tol tolerance for optimizaiton stop
#' @param prior indicates what penalty to use on the likelihood, if any
#' @return list of results
#' \item{V}{the estimated correlation matrix}
#' \item{status}{whether the correlation is global or local max}
#' \item{ll}{the loglikelihood value at different starting point}
#' \item{rho}{the estimated correlation at different stating point}
#' \item{ttime}{the total running time for each entry}
estimateV = function(data, Ulist, gridmult= sqrt(2),
                     grid = NULL,
                     normalizeU = TRUE,
                     usepointmass = TRUE, init_V = list(diag(ncol(data$Bhat))),
                     init_pi=NULL, init_rho=0, tol=1e-5, prior=c("nullbiased", "uniform"),
                     optmethod = c('optim', 'em', 'mrho', 'mV')){
  optmethod = match.arg(optmethod)
  Z = data$Bhat/data$Shat

  p = ncol(Z)
  V = diag(p)
  status = c()
  ll = list()
  rho = list()
  ttime = c()
  if(optmethod == 'mV'){
    result = mixture.MV.times(Z, Ulist, init_V = init_V, tol=tol, prior = prior)
    return(list(V = result$result$V, status = result$status,
                ll=result$loglik, ttime = result$time))
  }
  else{
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        # pairwise data
        Z.ij = Z[,c(i,j)]
        Ulist.ij = pairwise_cov(Ulist, i,j)

        if(optmethod == 'optim'){
          m.Z = mash_set_data(Z.ij, 1)
          if(normalizeU){Ulist.ij = mashr:::normalize_Ulist(Ulist.ij)}
          if(missing(grid)){
            grid.ij = mashr:::autoselect_grid(m.Z, gridmult)
            xUlist.ij = mashr:::expand_cov(Ulist.ij, grid.ij, usepointmass)
          }else{
            xUlist.ij = mashr:::expand_cov(Ulist.ij, grid, usepointmass)
          }

          result = optimize_pi_rho_times(Z.ij, xUlist.ij, init_rho=init_rho, prior=prior, tol=tol)
        }else if(optmethod == 'em'){
          m.Z = mash_set_data(Z.ij, 1)
          if(normalizeU){Ulist.ij = mashr:::normalize_Ulist(Ulist.ij)}
          if(missing(grid)){
            grid.ij = mashr:::autoselect_grid(m.Z, gridmult)
            xUlist.ij = mashr:::expand_cov(Ulist.ij, grid.ij, usepointmass)
          }else{
            xUlist.ij = mashr:::expand_cov(Ulist.ij, grid, usepointmass)
          }

          result = mixture.EM.times(Z.ij, xUlist.ij, init_pi=init_pi, init_rho=init_rho, prior=prior)
        }else if(optmethod == 'mrho'){
          result = mixture.M.rho.times(Z.ij, Ulist.ij, init_rho = init_rho, prior=prior, tol=tol)
        }
        V[i,j] = result$result$rho
        V[j,i] = V[i,j]
        status = c(status, result$status)
        ll = c(ll, list(result$loglik))
        rho = c(rho, list(result$rho))
        ttime = c(ttime, sum(result$time))
      }
    }
    return(list(V = V, status = status, ll = ll, rho=rho, ttime = ttime))
  }
}

pairwise_cov = function(Ulist, i,j){
  lapply(Ulist, function(U) U[c(i,j), c(i,j)])
}
