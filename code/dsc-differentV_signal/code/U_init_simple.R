B1 = matrix(0,n/2,p)
B2 = matrix(0,n/2,p)
B2[,1] = rnorm(n/2)
B = rbind(B1, B2)
Vtrue = array(0, dim=c(p,p,n))
Ehat = matrix(0,n,p)
for(i in 1:n){
  Vtrue[,,i] = clusterGeneration::rcorrmatrix(p)
  Ehat[i,] = mvtnorm::rmvnorm(1, sigma = Vtrue[,,i])
}
Bhat = B + Ehat
data = list(B = B, Bhat=Bhat, Shat = 1)
