lambda <- function(x,y,t) {
  one <- (t^2)*b1*dmvnorm(matrix(c(x,y),ncol=2), mu1, s1)
  two <- (1-t)*b2*dmvnorm(matrix(c(x,y),ncol=2), mu2, s2)
  three <- t*(1-t)*b3*dmvnorm(matrix(c(x,y),ncol=2), mu3, s3)
  lambda_0*exp(b0+one+two+three)
}

lambda_e <- function(x,y){
  lambda(x,y,point)
}

lambda_t <- function(x,y,t) {
  one <- 2*t*b1*dmvnorm(c(x,y), mu1, s1)
  two <- -b2*dmvnorm(c(x,y), mu2, s2)
  three <- (1-2*t)*b3*dmvnorm(c(x,y), mu3, s3)
  lambda(x,y,t)*(one+two+three)
}

derx1 <- function(x,y,t){
  dmvnorm(c(x,y), mu1, s1)*(-(x-mu1[1])*si1[1,1]-(y-mu1[2])*si1[1,2])
}

derx2 <- function(x,y,t){
  dmvnorm(c(x,y), mu2, s2)*(-(x-mu2[1])*si2[1,1]-(y-mu2[2])*si2[1,2])
}

derx3 <- function(x,y,t){
  dmvnorm(c(x,y), mu3, s3)*(-(x-mu3[1])*si3[1,1]-(y-mu3[2])*si3[1,2])
}

lambda_x <- function(x,y,t) {
  one <- t*t*b1*derx1(x,y,t)
  two <- (1-t)*b2*derx2(x,y,t)
  three <- t*(1-t)*b3*derx3(x,y,t)
  lambda(x,y,t)*(one+two+three)
}

dery1 <- function(x,y,t){
  dmvnorm(c(x,y), mu1, s1)*(-(x-mu1[1])*si1[1,2]-(y-mu1[2])*si1[2,2])
}

dery2 <- function(x,y,t){
  dmvnorm(c(x,y), mu2, s2)*(-(x-mu2[1])*si2[1,2]-(y-mu2[2])*si2[2,2])
}

dery3 <- function(x,y,t){
  dmvnorm(c(x,y), mu3, s3)*(-(x-mu3[1])*si3[1,2]-(y-mu3[2])*si3[2,2])
}

lambda_y <- function(x,y,t) {
  one <- t*t*b1*dery1(x,y,t)
  two <- (1-t)*b2*dery2(x,y,t)
  three <- t*(1-t)*b3*dery3(x,y,t)
  lambda(x,y,t)*(one+two+three)
}

Norma <- function(x,y) {sqrt(x**2+y**2)}

true.values <- function(dp,time,per){
  dp.new <- dp[(dp[,1]>=0)&(dp[,1]<=1)&(dp[,2]>=0)&(dp[,2]<=1),]
  inten <- mapply(lambda,dp.new[,1],dp.new[,2],time)
  dert <- mapply(lambda_t,dp.new[,1],dp.new[,2],time)
  derx <- mapply(lambda_x,dp.new[,1],dp.new[,2],time)
  dery <- mapply(lambda_y,dp.new[,1],dp.new[,2],time)
  norm_grad <- mapply(Norma,derx,dery)
  
  vel <- abs(dert)/norm_grad
  a1 <- data.frame(x=dp.new[,1],y=dp.new[,2])
  
  h.x <- unique(a1$x)[2]-unique(a1$x)[1]
  h.y <- unique(a1$y)[2]-unique(a1$y)[1]
  
  a1[,"inten"] <- inten
  a1[,"difx_pre"] <- derx
  a1[,"dify_pre"] <- dery
  a1[,"dt"] <- dert
  a1[,"vel"] <- vel
  a1["norm"] <- norm_grad
  a1["difxn"] <- (a1["difx_pre"]/a1$norm)*h.x
  a1["difyn"] <- (a1["dify_pre"]/a1$norm)*h.y
  a1["difxn2"] <- sign(a1$dt)*a1["difxn"]
  a1["difyn2"] <- sign(a1$dt)*a1["difyn"]
  
  a1[,"velt"] <- a1$vel
  sup <- round(quantile(a1$vel,per,names=FALSE,na.rm=TRUE),1)
  a1[!is.na(a1$vel) & a1$vel>sup,"velt"] <- sup*0.99
  return(a1)
}