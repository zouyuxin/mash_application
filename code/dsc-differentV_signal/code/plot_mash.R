pdf(debug_plots)
# corrplot::corrplot(V$mash.model$fitted_g$Ulist[['ED_PCA_1']], method='color', cl.lim=c(-1,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45,
#                    col=colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF", "#E0F3F8","#91BFDB","#4575B4")))(128))
# corrplot::corrplot(V$mash.model$fitted_g$Ulist[['ED_PCA_2']], method='color', cl.lim=c(-1,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45,
#                    col=colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF", "#E0F3F8","#91BFDB","#4575B4")))(128))
barplot(mashr::get_estimated_pi(V$mash.model), las=2, cex.names = 0.7)
if(all(!is.na(mashr::get_pairwise_sharing(V$mash.model)))){
  corrplot::corrplot(mashr::get_pairwise_sharing(V$mash.model), method='color', cl.lim=c(-1,1), type='upper', addCoef.col = "black", tl.col="black", tl.srt=45,
                     col=colorRampPalette(rev(c("#D73027","#FC8D59","#FEE090","#FFFFBF", "#E0F3F8","#91BFDB","#4575B4")))(128))
}
dev.off()
