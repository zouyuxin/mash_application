library(mashr)
# mash cov
m.data = mash_set_data(Bhat = data$Bhat, Shat = data$Shat)
m.1by1 = mash_1by1(m.data)

# strong = get_significant_results(m.1by1)
# U.pca = cov_pca(m.data, min(5,p), subset = strong)
# U.ed = cov_ed(m.data, U.pca, subset = strong)
U.c = cov_canonical(m.data)
# Ulist = c(U.c, U.ed)
Ulist = c(U.c)
