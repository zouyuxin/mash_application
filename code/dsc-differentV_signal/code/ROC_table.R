ROC.table = function(data, model){
  sign.test = data*model$result$PosteriorMean
  thresh.seq = seq(0, 1, by=0.005)[-1]
  m.seq = matrix(0,length(thresh.seq), 2)
  colnames(m.seq) = c('TPR', 'FPR')
  for(t in 1:length(thresh.seq)){
    m.seq[t,] = c(sum(sign.test>0 & model$result$lfsr <= thresh.seq[t])/sum(data!=0),
                  sum(data==0 & model$result$lfsr <=thresh.seq[t])/sum(data==0))
  }
  return(m.seq)
}
