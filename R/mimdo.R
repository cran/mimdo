#' @return The output returns a data frame of the complete imputed data. This means that the missing values of the original incomplete dataset have been imputed. If the function does not return a value, this means that the covariance matrix is not invertible and is exactly singular.
#' @export
#' @importFrom stats optim
#' @importFrom stats cov
#' @param incomplete_data A data frame with missing values.
#' @param inverse If TRUE, the inverse covariance matrix will be used for distance calculation. If the covariance matrix is non-invertible, use inverse = FALSE.
#' @param iterations Number of iterations. It can be adjusted to avoid long running time.
#' @title Multivariate Imputation by Mahalanobis Distance Optimization
#' @description
#' Imputes missing values of an incomplete data matrix by minimizing the Mahalanobis distance of each sample from the overall mean. By utilizing Mahalanobis distance, this imputation method is preferable to be used on datasets with highly correlated variables.
#' @author Geovert John D. Labita
#' @references Labita, GJ.D. and Tubo, B.F. (2024). Missing data imputation via optimization approach: An application to K-means clustering of extreme temperature. Reliability: Theory and Applications, 2(78), 115-123. DOI: https://doi.org/10.24412/1932-2321-2024-278-115-123
#' @references Bertsimas, D., Pawlowski, C., and Zhou, Y.D. (2018). From predictive methods to missing data imputation: An optimization approach. Journal of Machine Learning Research, 18(196), 1-39.
#' @examples
#' incomplete_data<-as.data.frame(matrix(c(5.1,NA,4.7,NA,3.0,3.2,1.4,1.4,NA,0.2,0.2,NA),nrow=3))
#' mimdo(incomplete_data, inverse=FALSE)
mimdo<-function(incomplete_data,inverse,iterations=30){

  #indices of missing data
  d2<-dim(incomplete_data)
  complete_data<-incomplete_data
  for (i in 1:d2[1]){
    for (j in 1:d2[2]){
      if (is.na(complete_data[i,j])==TRUE){
        complete_data[i,j]=mean(complete_data[,j], na.rm=TRUE)}}}
  datum1<-incomplete_data
  indices<-numeric()
  e2<-numeric()
  for (i in 1:d2[1]){
    for (j in 1:d2[2]){
      if (is.na(incomplete_data[i,j])==TRUE){
        indices<-append(indices,i)
        e2<-append(e2,j)}}}
  special<-unique(indices)
  d3<-length(special)
  e1<-numeric()
  imputed<-incomplete_data[rowSums(is.na(incomplete_data))>0,]
  for (i in 1:d3){
    for (j in 1:d2[2]){
      if (is.na(imputed[i,j])==TRUE){
        e1<-append(e1,i)}}}

  #imputation
  resultant<-function(input,input2){
    initial<-as.matrix(colMeans(input))
    input2=inverse
    if (input2==TRUE){
      covar<-solve(cov(input))}
    if (input2==FALSE){
      covar<-cov(input)}
    minus<-matrix(nrow=d2[1], ncol=d2[2])
    plus<-matrix(nrow=d2[1], ncol=d2[2])
    for (i in 1:d2[1]){
      for (j in 1:d2[2]){
        minus[i,j]<-datum1[i,j]-initial[j,1]
        plus[i,j]<-input[i,j]-initial[j,1]}}
    for (i in 1:d2[1]){
      for (j in 1:d2[2]){
        if (is.na(datum1[i,j])==TRUE){
          minus[i,j]=0}}}
    product1<-matrix(nrow=d2[1], ncol=d2[2])
    product2<-matrix(nrow=d2[1], ncol=d2[2])
    for (c in 1:d2[1]){
      for (b in 1:d2[2]){
        u1=minus[c,]*covar[,b]
        product1[c,b]=sum(u1)
        u2=plus[c,]*covar[,b]
        product2[c,b]=sum(u2)}}
    mahalanobis<-numeric()
    connect=matrix(nrow=0, ncol=d2[2])
    for (j in 1:d2[1]){
      v1=product1[j,]*minus[j,]
      connect<-rbind(connect,v1)
      v2=product2[j,]*plus[j,]
      mahalanobis<-append(mahalanobis,sum(v2))}
    connection<-matrix(nrow=d3, ncol=d2[2])
    c0<-matrix(nrow=d3, ncol=d2[2])
    for (i in 1:d3){
      connection[i,]=connect[special[i],]
      c0[i,]=minus[special[i],]}
    dim<-length(indices)
    wire<-matrix(nrow=0, ncol=d2[2])
    c2<-matrix(nrow=0, ncol=d2[2])
    for (j in 1:dim){
      wire=rbind(wire,connection[e1[j],])
      c2<-rbind(c2,c0[e1[j],])}
    c1<-rowSums(wire)
    outcome1<-numeric()
    for (i in 1:dim){
      profit<-function(x){
        c3<-crossprod(c2[i,],covar[e2[i],])
        c4<-initial[e2[i],1]
        return(c1[i]+(x-c4)*(covar[e2[i],e2[i]]*(x-c4)+c3))
      }
      c5<-optim(1,profit,lower=0,method="L-BFGS-B")
      outcome1<-append(outcome1,c5$par)}
    return(outcome1)}

  #iteration
  stop<-numeric()
  stop<-append(stop,1)
  output<-resultant(complete_data)
  out<-length(output)
  outputs<-matrix(nrow=out,ncol=0)
  outputs<-cbind(outputs,output)
  primary<-incomplete_data
  j=1
  for (a in 1:d2[1]){
    for (b in 1:d2[2]){
      if (is.na(primary[a,b])==TRUE){
        primary[a,b]=output[j]
        j=j+1}}}
  while (length(stop)!=iterations){
    stop<-append(stop,1)
    output<-resultant(primary)
    outputs<-cbind(outputs,output)
    primary<-incomplete_data
    j=1
    for (a in 1:d2[1]){
      for (b in 1:d2[2]){
        if (is.na(primary[a,b])==TRUE){
          primary[a,b]=output[j]
          j=j+1}}}}
  goal<-rowMeans(outputs)
  primary<-incomplete_data
  j=1
  for (a in 1:d2[1]){
    for (b in 1:d2[2]){
      if (is.na(primary[a,b])==TRUE){
        primary[a,b]=goal[j]
        j=j+1}}}
  return(primary)}
