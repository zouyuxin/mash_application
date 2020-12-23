B = matrix(0,n,p)
Vtrue = array(0, dim=c(p,p,n))
Bhat = matrix(0,n,p)
for(i in 1:n){
  Vtrue[,,i] = clusterGeneration::rcorrmatrix(p)
  Bhat[i,] = mvtnorm::rmvnorm(1, sigma = Vtrue[,,i])
}
data = list(B = B, Bhat=Bhat, Shat = 1)
