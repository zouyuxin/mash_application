generate_data = function(n, p, V, Utrue, pi=NULL){
  if (is.null(pi)) {
    pi = rep(1, length(Utrue)) # default to uniform distribution
  }
  assertthat::are_equal(length(pi), length(Utrue))

  for (j in 1:length(Utrue)) {
    assertthat::are_equal(dim(Utrue[j]), c(p, p))
  }

  pi <- pi / sum(pi) # normalize pi to sum to one
  which_U <- sample(1:length(pi), n, replace=TRUE, prob=pi)

  Beta = matrix(0, nrow=n, ncol=p)
  for(i in 1:n){
    Beta[i,] = MASS::mvrnorm(1, rep(0, p), Utrue[[which_U[i]]])
  }

  E = MASS::mvrnorm(n, rep(0, p), V)
  Bhat = Beta + E
  Shat = 1
  return(list(B = Beta, Bhat=Bhat, Shat = Shat))
}

V.simulation = function(n,p,times,Utrue, seed = 1, optmethod = c('mle', 'em1', 'em2')){
  set.seed(seed)
  ttime = c()
  mFE.max = numeric(times)
  mFE.max.nearPD = numeric(times)
  mFE.trun = numeric(times)
  m2E.max = numeric(times)
  m2E.max.nearPD = numeric(times)
  m2E.trun = numeric(times)
  mashlik.true = numeric(times)
  mashlik.Vhat = numeric(times)
  mashlik.Vhat.nearF = numeric(times)
  mashlik.Vhat.near2 = numeric(times)
  Vt = list()
  Vmax = list()
  Vtrun = list()
  pd = 0
  status.g = 0
  for(t in 1:times){
    Vtrue = clusterGeneration::rcorrmatrix(p)
    Vt = c(Vt, list(Vtrue))
    data = generate_data(n, p, Vtrue, Utrue)
    # mash cov
    m.data = mash_set_data(Bhat = data$Bhat, Shat = data$Shat)
    m.1by1 = mash_1by1(m.data)
    strong = get_significant_results(m.1by1)

    U.pca = cov_pca(m.data, 3, subset = strong)
    U.ed = cov_ed(m.data, U.pca, subset = strong)
    U.c = cov_canonical(m.data)

    Vhat.max <- estimateV(m.data, c(U.c, U.ed), init_rho = c(-0.5,0,0.5), tol=1e-4, optmethod = optmethod)
    ttime = c(ttime, sum(Vhat.max$ttime))

    Vmax = c(Vmax, list(Vhat.max$V))

    m.data.true = mash_set_data(Bhat = data$Bhat, Shat = data$Shat, V = Vtrue)
    m.model.true = mash(m.data.true, c(U.c,U.ed), verbose = FALSE)
    mashlik.true[t] = get_loglik(m.model.true)

    # check pd
    R <- tryCatch(chol(Vhat.max$V),error = function (e) FALSE)
    if (is.matrix(R)){
      pd = pd + 1
      m.data.Vhat = mash_set_data(Bhat = data$Bhat, Shat = data$Shat, V = Vhat.max$V)
      m.model.Vhat = mash(m.data.Vhat, c(U.c,U.ed), verbose = FALSE)
      mashlik.Vhat[t] = get_loglik(m.model.Vhat)
    }else{
      Vhat.near.F = as.matrix(Matrix::nearPD(Vhat.max$V, conv.norm.type = 'F', keepDiag = TRUE)$mat)
      Vhat.near.2 = as.matrix(Matrix::nearPD(Vhat.max$V, conv.norm.type = '2', keepDiag = TRUE)$mat)

      # Error
      mFE.max.nearPD[t] = norm(Vhat.near.F - Vtrue, 'F')
      m2E.max.nearPD[t] = norm(Vhat.near.2 - Vtrue, '2')

      # likelihood
      m.data.Vhat = mash_set_data(Bhat = data$Bhat, Shat = data$Shat, V = Vhat.near.F)
      m.model.Vhat = mash(m.data.Vhat, c(U.c,U.ed), verbose = FALSE)
      mashlik.Vhat.nearF[t] = get_loglik(m.model.Vhat)

      m.data.Vhat = mash_set_data(Bhat = data$Bhat, Shat = data$Shat, V = Vhat.near.2)
      m.model.Vhat = mash(m.data.Vhat, c(U.c,U.ed), verbose = FALSE)
      mashlik.Vhat.near2[t] = get_loglik(m.model.Vhat)
    }

    # check global
    status.g = status.g + sum(Vhat.max$status == 'global')
    # compute error
    mFE.max[t] = norm(Vhat.max$V - Vtrue, 'F')
    m2E.max[t] = norm(Vhat.max$V - Vtrue, '2')

    # truncated cor
    Vhat.tru = estimate_null_correlation(m.data)
    Vtrun = c(Vtrun, list(Vhat.tru))
    mFE.trun[t] = norm(Vhat.tru - Vtrue, 'F')
    m2E.trun[t] = norm(Vhat.tru - Vtrue, '2')
  }
  return(list(Time = ttime, mFE.max = mFE.max, mFE.trun = mFE.trun,
              m2E.max = m2E.max, m2E.trun = m2E.trun,
              mFE.max.nearPD = mFE.max.nearPD, m2E.max.nearPD = m2E.max.nearPD,
              Vt = Vt, Vmax = Vmax, Vtrun = Vtrun, pd = pd, status.g = status.g,
              mashloglik.true = mashlik.true, mashloglik.Vhat = mashlik.Vhat,
              mashloglik.Vhat.nearF = mashlik.Vhat.nearF, mashloglik.Vhat.near2 = mashlik.Vhat.near2))
}
