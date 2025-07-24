## This code has all the functions needed to estimate the derivative
## of a function using regression splines. 
## This code was created in on the 9th of January 2025.


## Load the required libraries
library(mgcv)

## The lines pf codes below defines the base functions for the paper###############
tpower <- function(x, t, p)
{
  ## Truncated p-th power function
  ## x is the values to be evaluated at knot t
  ## p is the degree of the plynomial
  
  ((x-t)^p)*(x>t)
}

### Bspline function
Bbase <- function (x, xl = min(x), xr = max(x), nseg = 10, bdeg = 3) 
{
  if (xl > min(x)) {
    xl = min(x)
    warning("Left boundary adjusted to min(x) = ", xl)
  }
  if (xr < max(x)) {
    xr = max(x)
    warning("Right boundary adjusted to max(x) = ", xr)
  }
  dx <- (xr - xl)/nseg
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1) * dx^bdeg)
  B <- (-1)^(bdeg + 1) * P %*% t(D)
  dim(B)
  # par(mfrow=c(1,2))
  # matplot(x,B, type = "l", lwd = 2)
  nb <- ncol(B)
  sk <- knots[(1:nb) + bdeg + 1]
  Mask <- outer(x, sk, "<")
  B <- B * Mask
  # dim(B1)
  # matplot(x,B1, type = "l", lwd = 2)
  
  att2 = list(x = x, xl = xl, xr = xr, nseg = nseg, bdeg = bdeg, 
              B=B, knots=knots)
  return(att2)
}
# ' @title B-spline basis creation function
## Bspline basis creation function on a grid.
bbase.grid <- function(x,dx,knots,bdeg)
{
  P <- outer(x,knots,tpower,bdeg)
  ## P is the basis matrix
  n <- dim(P)[2]

  D <- diff(diag(n),diff=bdeg+1)/(gamma(bdeg+1)*dx^bdeg)

  B <- (-1)^(bdeg+1)*P %*% t(D)
  
  return(B)
}


## This is the main function to estimate the derivative function
## using the penalized spline method. We use decomaposition approach to enhance comqputational efficiency.
pgams <- function(x=x,y=y,lambda=0.1, r=0, x.grid, nseg=35,pord=2,bdeg=3){
  n <- length(x)
  BS<- Bbase(x, nseg = nseg, bdeg = bdeg)
  BS.grid<- Bbase(x.grid, nseg = nseg, bdeg = bdeg)
  B.grid <-BS.grid$B
  dx <- (BS$xr- BS$xl)/nseg
  knots <- BS$knots
  B <- BS$B
  ####. Cholesky Decomposition
  R <- chol(t(B)%*%B)
  ### making my D matrix
  v=dim(BS$B)[2]
  Dm <- diff(diag(v), diff = pord)
  Mv <- solve(t(R))%*%t(Dm)%*%Dm%*%solve(R)
  ## Performing svd to find the Gamma and T inverse 
  ssvd <- svd(Mv)
  Tinves <- solve(R)%*%ssvd$u
  Gamma <- diag(ssvd$d)
  #gama<-mean(diag(Gamma))
  ### Finding reparametrised parameters
  hK <- solve(diag(1, nrow = v)+lambda*Gamma)%*%(t(B%*%Tinves))
  Thetahat <- hK%*%y
  A <- B%*%Tinves%*%hK
  f.hat <- A%*%y
  
  #theta <-1/(1+lambda*gama)*diag(1, nrow = v)%*%t(B%*%Tinves)%*%y
  ## Extract the alpha hat(original estimates)
  alpha.hat <-  Tinves %*% Thetahat
  fg.hat <- B.grid%*%alpha.hat
  ##### Fitting the derivative function
  if(r==0){
    Dr <- diag(v)
  }else{
    Dr <- diff(diag(v), differences = r)
    
  }
  Bqrx<- Bbase(x, nseg = nseg, bdeg = bdeg-r)
  Bqrx.grid<- Bbase(x.grid, nseg = nseg, bdeg = bdeg-r)
  
  Atilde<- Bqrx$B%*%Dr%*%Tinves
  Bq = (Bqrx$B%*%Dr)/(dx^r)
  Bq.grid = (Bqrx.grid$B%*%Dr)/(dx^r)
  K <- Bq%*%Tinves%*%hK
  fr.hat <- Bq%*%alpha.hat
  frg.hat<- Bq.grid%*%alpha.hat
  return(list(x.grid=x.grid,BS=BS, f.hat=f.hat,fg.hat=fg.hat,
              fr.hat=fr.hat,frg.hat, K=K, M=hK, Atilde=Atilde,
              A=A, lambda=lambda))
  
}


##### Naive derivative estimation using GCV to choose the smoothing parameter ##################

### Generalised cross validation criterion
## GCV smoothing parameter selector for Naive estimator
gcvlambda <- function(lambda = 0,x=x, y=y,nseg=35, pord = 3,
                      bdeg=4){
  
  fit0 <- pgams(x=x,y=y, lambda = lambda, r=0, x.grid = x, nseg = nseg, pord = pord,
        bdeg = bdeg)
  n <- length(x)
  RSS_lambda<- as.numeric((t(y-fit0$f.hat)%*%(y-fit0$f.hat)))
  EDF_lambda <- sum(diag(fit0$A))
  
  logGCV_lambda <- log(RSS_lambda)-2*log((n/sqrt(n))-EDF_lambda/sqrt(n) )
  
  return(logGCV_lambda)
}


## Naive estimate of the mean and derivative function. This uses grid search
## to find the optimal smoothing parameter  
naive.est <- function(x=x, y=y,r=r,log.lam= seq(-2.5, 2, length.out=800)
                      , nseg = 35, bdeg = 4,pord = 2,
                      x.grid = x.grid){
  n=length(x)
  lm1<- 10^log.lam
  searh<- sapply(lm1, function(l){
    gcvlambda(lambda = l, x=x,y=y,
              nseg=nseg, pord = pord, bdeg = bdeg)
  })
  
  
  idx<- which.min(searh)
  fr.est <- pgams(x=x,y=y, lambda = lm1[idx], r=r, x.grid = x.grid, nseg = nseg, 
                  pord = pord,
                  bdeg = bdeg)
  sig.hat <- sqrt(sum ((y-fr.est$f.hat)^2)/(n-sum(diag(fr.est$A))))
  return(list(fr.est=fr.est, f.hat= fr.est$f.hat,fg.hat=fr.est$fg.hat,
              fr.hat= fr.est$fr.hat,frg.hat= fr.est$frg.hat,sig.hat=sig.hat,
              lambda=lm1[idx]))
}

## Naive estimate of the mean and derivative function. This uses optimization
## to find the optimal smoothing parameter
naive.est.opt <- function(x=x, y=y,r=r
                      , nseg = 35, bdeg = 4,pord = 2,
                      x.grid = x.grid){
  n <- length(x)
  initial.lambda <- 1
  searh<- optim(par = initial.lambda,fn=gcvlambda ,x=x,y=y,
        nseg=nseg, pord = pord,
        bdeg=bdeg,
        method = "L-BFGS-B",lower=0,upper=Inf)
  
  idx<- searh$par
  fr.est <- pgams(x=x,y=y, lambda = idx, r=r, x.grid = x.grid, nseg = nseg, 
                  pord = pord,
                  bdeg = bdeg)
  edf<- (n-sum(diag(fr.est$A)))
  rss <- sum((y-fr.est$f.hat)^2)
  sig.hat <- sqrt((rss/edf) )
  BS <- fr.est$BS
  return(list(BS=BS, fr.est=fr.est, f.hat= fr.est$f.hat,fg.hat=fr.est$fg.hat,
              fr.hat= fr.est$fr.hat,frg.hat= fr.est$frg.hat,sig.hat=sig.hat,
              lambda=idx,edf=edf, tr=sum(diag(fr.est$A))))
}

########## End of naive estimation functions ######################

   






##### The proposed method for estimating the derivative functions ####################



#################################### The simplest estimate of fr ############################
simple.est <- function(lambda, r, x,x.grid, f,bdeg, nseg,pord){
  BS<- Bbase(x, nseg = nseg, bdeg = bdeg)
  dx <- (BS$xr- BS$xl)/nseg
  bdeg <- BS$bdeg
  knots <- BS$knots
  B <- BS$B
  v <- dim(B)[2]
  Dm <- diff(diag(v), differences = pord)
  Pm <- t(Dm) %*% Dm
  Dr <- diff(diag(v), differences =r)
  B.qr<- Bbase(x, nseg = nseg, bdeg = bdeg-r)
  Btil <-as.matrix((B.qr$B%*%Dr)/(dx^r))
  
  H <- Btil %*% solve(t(B)%*%B + lambda*Pm)%*%t(B)
  B.qr.grid<- Bbase(x.grid, nseg = nseg, bdeg = bdeg-r)
  Btil.grid <-as.matrix((B.qr.grid$B%*%Dr)/(dx^r))
  H.grid<- Btil.grid %*% solve(t(B)%*%B + lambda*Pm)%*%t(B)
  # H2<- Btil %*% solve(t(B)%*%B + lambda2*Pm)%*%t(B) 
  H.f <- as.vector(H%*%f)
  H.f.grid <- as.vector(H.grid%*%f)
  return(list(x=x,fr.hat = H.f,fr.hat.grid = H.f.grid))
}

### The objective function
mise.lambda <- function(lambda=0.1, x,y, r=1,sig=0.1,nseg=35, pord=2, 
                    bdeg=35, f=f, fr=NULL)
{
  BS<- Bbase(x, nseg = nseg, bdeg = bdeg)
  dx <- (BS$xr- BS$xl)/nseg
  bdeg <- BS$bdeg
  knots <- BS$knots
  B <- BS$B
  v <- dim(B)[2]
  Dm <- diff(diag(v), differences = pord)
  Pm <- t(Dm) %*% Dm
  Dr <- diff(diag(v), differences =r)
  B.qr<- Bbase(x, nseg = nseg, bdeg = bdeg-r)
  Btil <-as.matrix((B.qr$B%*%Dr)/(dx^r))
  A.prime.lam <- Btil %*% solve(t(B)%*%B + lambda*Pm)%*%t(B) 
  A.prime.lam.f <- as.vector(A.prime.lam%*%f)
  
  var <- ((sig^2)*sum(A.prime.lam^2))/length(x)
  if(is.null(fr)){
    fr.est <- A.prime.lam%*%y
    bias <- (A.prime.lam.f- fr.est)
  }else{
    bias <- (A.prime.lam.f- fr)
  }
  
  # bias <- (A.prime.lam.f- fr)
  sq.bias<- sum((bias^2))/length(x)
  
  return( list(mise=(var+sq.bias), var= var, sq.bias=sq.bias,
               H=A.prime.lam) )
}


## This is the objective function for the optimization
## For optimization. This is the estimation of the MISE
mise.lambda.optim <- function(lambda=0.1, x,y, r=1,sig=0.1,nseg=35, pord=2, 
                        bdeg=35, f=f, fr=NULL)
{
  BS<- Bbase(x, nseg = nseg, bdeg = bdeg)
  dx <- (BS$xr- BS$xl)/nseg
  bdeg <- BS$bdeg
  knots <- BS$knots
  B <- BS$B
  v <- dim(B)[2]
  Dm <- diff(diag(v), differences = pord)
  Pm <- t(Dm) %*% Dm
  Dr <- diff(diag(v), differences =r)
  B.qr<- Bbase(x, nseg = nseg, bdeg = bdeg-r)
  Btil <-as.matrix((B.qr$B%*%Dr)/(dx^r))
  A.prime.lam <- Btil %*% solve(t(B)%*%B + lambda*Pm)%*%t(B) 
  A.prime.lam.f <- as.vector(A.prime.lam%*%f)
  
  var <- ((sig^2)*sum(A.prime.lam^2))/length(x)
  if(is.null(fr)){
    fr.est <- A.prime.lam%*%y
    bias <- (A.prime.lam.f- fr.est)
  }else{
    bias <- (A.prime.lam.f- fr)
  }
  
  # bias <- (A.prime.lam.f- fr)
  sq.bias<- sum(bias^2)/length(x)
  
  return(mise=(var+sq.bias) )
}


## This is another variant of the objective function for the optimization
mise.lambda.optimVB <- function(lambda=0.1, x,y, r=1,sig=0.1,nseg=35, pord=2, 
                              bdeg=35, f=f, fr=NULL, p=2)
{
  n<- length(x)
  BS<- Bbase(x, nseg = nseg, bdeg = bdeg)
  dx <- (BS$xr- BS$xl)/nseg
  bdeg <- BS$bdeg
  knots <- BS$knots
  B <- BS$B
  v <- dim(B)[2]
  Dm <- diff(diag(v), differences = pord)
  Pm <- t(Dm) %*% Dm
  Dr <- diff(diag(v), differences =r)
  B.qr<- Bbase(x, nseg = nseg, bdeg = bdeg-r)
  Btil <-as.matrix((B.qr$B%*%Dr)/(dx^r))
  dim(Btil)
  A.prime.lam <- Btil %*% solve((t(B)%*%B + lambda*Pm))%*%t(B) 
  
  #A.prime.lamp <- Btil[p,] %*% solve((t(B)%*%B + lambda*Pm))%*%t(B) 
  
#   dim(A.prime.lam)
#   dim(A.prime.lamp)
#   
# all.equal(A.prime.lam[p,] ,as.vector(A.prime.lamp))

A.prime.lam.f <- as.vector(A.prime.lam%*%f)
  # dim(A.prime.lamp)
  b <- sum(A.prime.lam[p,]^2)
  b1 <- sum(A.prime.lam^2)
var <- (((sig^2)/n)*b1)
varp <- (((sig^2)/n)*b)
  if(is.null(fr)){
    fr.est <- A.prime.lam%*%y
    bias <- (A.prime.lam.f- fr.est)
  }else{
    bias <- (A.prime.lam.f- fr)
  }
  
  # bias <- (A.prime.lam.f- fr)
  sq.bias<- sum(bias^2)/n
  
  return(list( mise=(var+sq.bias), var=var,varp=varp, bias=mean(bias),
               HH = (b1^2), hh=(b^2) ) )
}


  

## This is one-step plug-in estimate function. We optimize the mise.lambda.optim function 
## plug-in estimate function
plugin.est<-function(x=x, y=y, r=1,nseg=35, pord = 3,
                        bdeg=4, x.grid=x,fr=NULL){
  ### make the naive estimate
  naive.fr<- naive.est.opt(x=x, y=y, r=r,nseg=nseg, pord =pord ,
                         bdeg=bdeg, x.grid=x.grid)
  initial.lambda <- 1
  ##### make the plug-in
  sig.hat <- (naive.fr$sig.hat)
if(is.null(fr)){
  est.plugin <- optim(par = initial.lambda,fn=mise.lambda.optim ,x=x,y=y,
                      r=r,sig=sig.hat,
                      nseg=nseg, pord = pord,
                      bdeg=bdeg, f=naive.fr$f.hat
                      ,fr=NULL,
                      method = "L-BFGS-B",lower=0,upper=Inf)
}else{
  est.plugin <- optim(par = initial.lambda,fn=mise.lambda.optim ,x=x,y=y,
                       r=r,sig=sig.hat,
                       nseg=nseg, pord = pord,
                       bdeg=bdeg, f=naive.fr$f.hat,
                       fr=naive.fr$fr.hat,
                       method = "L-BFGS-B",lower=0,upper=Inf)
}
  lambda_plugin <- (est.plugin$par)
  plugin.fr<- pgams(lambda = lambda_plugin,x=x, y=y, r=r,nseg=nseg, 
                    pord =pord,bdeg=bdeg,x.grid = x.grid)
  
  return(list(x.grid=x.grid,f.hat=plugin.fr$f.hat,fg.hat=plugin.fr$fg.hat,
              fr.hat=plugin.fr$fr.hat,frg.hat=plugin.fr$frg.hat,
              lambda=lambda_plugin,
              K=plugin.fr$K, M=plugin.fr$M, Atilde=plugin.fr$Atilde,
              A=plugin.fr$A, sig.hat=sig.hat))
}



### This is the iterative re-substitution method function.  This is our main contribution
## to the paper. This is the iterative re-substitution method.
resub.est<- function(x=x, y=y, r=r, x.grid = x.grid,nseg=nseg, pord = pord,
                     bdeg=bdeg,tol=1e-10, ITs=10){
  
  keep <- matrix( NA,ITs,length(x) )
  keep.l <- numeric()
  plugin.fit <- plugin.est(x=x, y=y, r=r, nseg=nseg, 
                             pord = pord,bdeg=bdeg, x.grid = x.grid,fr="Y")
  
  fr.est <- plugin.fit 
  sig.hat <- plugin.fit$sig.hat
  f.hat <- fr.est$f.hat
  keep[1,] <- fr.est$fr.hat
  initial.lambda<- plugin.fit$lambda
  keep.l[1] <- initial.lambda
  for (i in 2:ITs)
  {
    est.resub <- optim(par = initial.lambda,fn=mise.lambda.optim,x=x,y=y,
                       r=r,sig=sig.hat,
                       nseg=nseg, pord = pord,
                       bdeg=bdeg, f=f.hat,
                       fr=fr.est$fr.hat,
                       method = "L-BFGS-B",lower=0,upper=Inf)
    
    lambda_resub <- (est.resub$par)
    update<- pgams(lambda = lambda_resub,x=x, y=y, r=r,nseg=nseg, 
                   pord =pord,
                   bdeg=bdeg,x.grid = x.grid)
    
    fr.est<- update
    f.hat <- fr.est$f.hat
    keep[i,] <- fr.est$fr.hat
    dif<- mean(abs(keep[i,]-keep[(i-1),]))
    keep.l[i]<- lambda_resub
    diff.l <- abs(keep.l[i]-keep.l[i-1])
    print(diff.l)
    print("*****")
    print(lambda_resub)
    print("*****")
    if(diff.l<=tol){
      break
    }
  }
  
  resub.lambda <- fr.est$lambda
  print(resub.lambda)
  return(list(x.grid=x.grid, fr.hat= fr.est$fr.hat, lambda=resub.lambda) )
  
}


## Function to make a confidence interval for the derivative function
#' @param deriv.est The estimated derivative function
#' @param r The order of the derivative
#' @param BS The basis functions object containing the basis matrix and other parameters
#' @param pord The polynomial order of the spline basis
#' @param sig.est The estimated standard deviation of the noise
#' @param lambda.est The estimated smoothing parameter
#' @param x The grid points where the derivative is estimated
#' @return A list containing the estimated derivative, lower and upper confidence intervals
pointwiseConfInt <- function(deriv.est,r,BS,pord, sig.est, lambda.est, x){

  B<- BS$B
  bdeg<- BS$bdeg
  knots<- BS$knots
  nseg <- BS$nseg
  dx <- (BS$xr- BS$xl)/nseg

  v <- dim(B)[2]
  Dm <- diff(diag(v), diff = pord)
  Pm <- t(Dm) %*% Dm
  print("The basis matrix B is:")
  print(head(B))
  print("********************")
  inv.part<- solve( (t(B)%*%B)+lambda.est*Pm)
  print("The inverse computation was successful!!!!.")
  # B.grid <- bbase.grid(x.grid,dx,knots,bdeg)
  Dr <- diff(diag(v), diff=r)
  # B.qr <- Bbase(x=x, nseg = nseg, bdeg = bdeg-r)$B
  B.qr <- bbase.grid(x= x, dx =dx, knots =knots[(r+1):(length(knots)-r)], bdeg= bdeg-r)
  est.var <- (sig.est^2)*(B.qr%*%Dr)%*%inv.part %*%t(B)%*%t((B.qr%*%Dr)%*%inv.part %*%t(B)) /(dx^(2*r) )
  
  est.st.dev<- sqrt(diag(est.var))
  
  pointwiseCI.lower <- deriv.est-2*est.st.dev
  pointwiseCI.upper <- deriv.est+2*est.st.dev

  return(list(CI.lower=pointwiseCI.lower,CI.upper=pointwiseCI.upper))
}

######### End of the proposed method functions ######################


############### Making the oracle estimator #######################
## This is the oracle estimator loss  function.
oracle.loss <- function(lambda=0.2, x, y, r, fr.grid, nseg=35, pord=2, bdeg=5, x.grid){
  fit <- pgams(lambda = lambda, x = x, y = y, r = r, x.grid = x.grid, nseg = nseg, pord = pord, bdeg = bdeg)
  loss.fun <- mean((fit$fr.hat - fr.grid)^2)
  return(loss.fun)
}

oracle.est <- function(initial.lambda = 0.03, x, y, r, fr.grid, nseg = 35, pord = 2, bdeg = 5, x.grid){
  est.oracle <- optim(par = initial.lambda, fn = oracle.loss, x = x, y = y, r = r, nseg = nseg, pord = pord, bdeg = bdeg, fr.grid = fr.grid, x.grid = x.grid, method = "L-BFGS-B", lower = 0, upper = Inf)
  lambda.oracle <- est.oracle$par
  model.oracle <- pgams(lambda = lambda.oracle, x = x, y = y, r = r, x.grid = x.grid, nseg = nseg, pord = pord, bdeg = bdeg)
  return(list(x.grid = x.grid, fr.hat= model.oracle$fr.hat, lambda = lambda.oracle))
}


### This is the oracle estimator loss  function.###########




########## Other functions #####################


## This function estimates the MISE for a given lambda and sigma
sig.effect<- function(sig, lams, x,r, f,fr,nseg = 35,pord = 2,
                      bdeg=4,log.lam= seq(-2.5, 1, length.out=800)){
  n<- length(x)
  y <- f + sig*rnorm(n)
  ## estimate of the mean function for a given lambda by GCV
  n.est <- naive.est(x=x, y=y,r=r,log.lam= log.lam
                     , nseg = nseg, bdeg = bdeg,pord = pord,
                     x.grid = x)
  ## This estimate the mise given fr.hat anf f.hat
  mise.est_ideal <- sapply(lams, function(l){ 
    mise.lambda(lambda = l, x=x,y=y, r=r, sig=n.est$sig.hat, nseg = nseg, 
                pord = pord,
                bdeg=bdeg, f=n.est$f.hat,fr=n.est$fr.hat)
  }
  )
  
  ## This estimate the mise given only f.hat 
  mise.est_prop <- sapply(lams, function(l){ 
    mise.lambda(lambda = l, x=x,y=y, r=r, sig=n.est$sig.hat, nseg = nseg, 
                pord = pord,
                bdeg=bdeg, f=n.est$f.hat, fr=NULL)
  }
  )
  
  mise.true <- sapply(lams, function(l){ 
    mise.lambda(lambda = l, x=x,y=y, r=r, sig=sig, nseg = nseg, pord = pord,
                bdeg=bdeg, f=f, fr=fr)
  }
  )
  return(list(mise.true=mise.true,mise.est_prop=mise.est_prop,
              mise.est_ideal=mise.est_ideal))
}






sig.lambda.fun <- function(sig,f,fr,x, x.grid, r, nseg, pord, bdeg 
                    ,tol=1e-10, ITs=10) {
  n<- length(x)
  y <- f + sig * rnorm(n)
  
  resub<- resub.est(x=x, y=y, r=r, x.grid = x.grid,nseg=nseg, pord = pord,
                      bdeg=bdeg,tol=tol, ITs=ITs)
  naive.fit <- naive.est.opt(x=x, y=y, r=r,nseg=nseg, pord = pord,
                         bdeg=bdeg, x.grid=x.grid)
  
  plugin.fitY <-  plugin.estfr(x=x, y=y, r=r,nseg=nseg, pord = pord,
                              bdeg=bdeg, x.grid=x.grid,fr="Y")
  plugin.fitN <-  plugin.estfr(x=x, y=y, r=r,nseg=nseg, pord = pord,
                             bdeg=bdeg, x.grid=x.grid,fr=NULL)
  
  initial.lambda <- plugin.fitY$lambda
  true.lambda <- optim(par = initial.lambda,fn=mise.lambda.optim ,x=x,y=y,
        r=r,sig=sig,
        nseg=nseg, pord = pord,
        bdeg=bdeg, f=f,
        fr=fr,
        method = "L-BFGS-B",lower=0,upper=Inf)
 
  
  return(list(true=true.lambda$par,naive=naive.fit$lambda,
              plugin.org=plugin.fitY$lambda,
              plugin.new=plugin.fitN$lambda, resub=resub$lambda ))
}


 







Jst<- function(x=x, y=y, r=r, x.grid = x.grid,nseg=nseg, pord = pord,
               bdeg=bdeg,tol=1e-10, ITs=10){
  
  keep <- matrix( NA,ITs,length(x) )
  keep.l <- numeric()
  plugin.fit <- plugin.estfr(x=x, y=y, r=r, nseg=nseg, 
                             pord = pord,bdeg=bdeg, x.grid = x.grid,fr="Y")
  
  fr.est <- plugin.fit 
  sig.hat <- plugin.fit$sig.hat
  f.hat <- fr.est$f.hat
  keep[1,] <- fr.est$fr.hat
  initial.lambda<- plugin.fit$lambda
  keep.l[1] <- initial.lambda
  Js<- numeric()
  Js[1]<- mise.lambda.optim(lambda =initial.lambda,x=x,y=y,
                            r=r,sig=sig.hat,
                            nseg=nseg, pord = pord,
                            bdeg=bdeg, f=f.hat,
                            fr=keep[1,] )
  for (i in 2:ITs)
  {
    est.resub <- optim(par = initial.lambda,fn=mise.lambda.optim,x=x,y=y,
                       r=r,sig=sig.hat,
                       nseg=nseg, pord = pord,
                       bdeg=bdeg, f=f.hat,
                       fr=fr.est$fr.hat,
                       method = "L-BFGS-B",lower=0,upper=Inf)
    
    lambda_resub <- (est.resub$par)
    update<- pgams(lambda = lambda_resub,x=x, y=y, r=r,nseg=nseg, 
                   pord =pord,
                   bdeg=bdeg,x.grid = x.grid)
    J<-  mise.lambda.optim(lambda =lambda_resub,x=x,y=y,
                           r=r,sig=sig.hat,
                           nseg=nseg, pord = pord,
                           bdeg=bdeg, f=f.hat,
                           fr=keep[i-1,] )
    Js[i]<- J
    fr.est<- update
    # f.hat <- fr.est$f.hat
    keep[i,] <- fr.est$fr.hat
    dif<- mean(abs(keep[i,]-keep[(i-1),]))
    keep.l[i]<- lambda_resub
    diff.l <- abs(keep.l[i]-keep.l[i-1])
    print(diff.l)
    print("*****")
    print(lambda_resub)
    print("*****")
    if(diff.l<=tol){
      break
    }
  }
  
  resub.lambda <- fr.est$lambda
  print(resub.lambda)
  return(list(x.grid=x.grid, fr.hat= fr.est$fr.hat, lambda=resub.lambda,
              lambdas =keep.l, Js= Js) )
  
}


## Function to compare the est mise and the true mise for a fixed lambda
samp <- function(range_start=100,range_end=1000,points_per_range=2 ){
  sample_size <- c()
  # Loop over each 1000-unit range and generate points
  for (i in seq(range_start, range_end - range_start, by = range_start)) {
    log_start <- log10(i)
    log_end <- log10(i + range_start)
    logn <- seq(log_start, log_end, length.out = points_per_range + 1)
    # Convert the log scale sequence to the normal scale and append to sample_size
    sample_size <- c(sample_size, round(10^(logn), 0))
  }
  # Ensure unique and sorted values
  sample_size <- sort(unique(sample_size))
  return(sample_size)
}


exp.mise <- function(lambda=0.2,r=1, bdeg=3, pord=2,sig=0.1,nsim=5,
                     log.n=seq(log10(1000), log10(10000),log10(1.11) )){
  library(dplyr)
  #ns <- samp(range_start=100,range_end=1000,points_per_range=l)
   v <- ( (2 * (2*pord))+ 1)
  # Ks <- ceiling(10*ns^((1/v)  ))
  # log.n<- seq(log10(1000), log10(10000),log10(1.11) )
  ns=ceiling(10^log.n)
  rate =1/v
  Ks <-ceiling( 10^( 1.6+rate*log.n) )
  
  keep.difmise <- data.frame(
    n = rep(NA, nsim*length(ns)),
    K = rep(NA, nsim*length(ns)),
    exp.vardiff= rep(NA, nsim*length(ns)),
    exp.vardiffp= rep(NA, nsim*length(ns)),
    exp.biasdiff= rep(NA, nsim*length(ns)),
    exp.sqdiff = rep(NA, nsim*length(ns)),
    sig.hats.diff = rep(NA, nsim*length(ns)),
    hh.sq = rep(NA, nsim*length(ns)),
    hh.sq.p = rep(NA, nsim*length(ns))
  )
  counter <- 1
  for(j in 1:length(ns)){
    #set.seed(200)
    x <- seq(0, 1, length.out=ns[j])
    f <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
    if(r==1){
      fr <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
    }else if(r==2){
      fr <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
    }
    p <- floor(ns[j]/2)
    for (i in 1:nsim) {
      
      y <- f + rnorm(n=ns[j],mean=0,sd=sig)
      
      t.emise <-  mise.lambda.optimVB(lambda=lambda, x=x,y=y, r=r,sig=sig,
                                    nseg=Ks[j], pord=pord, 
                                    bdeg=bdeg, f=f, fr=fr,p=p)
      
      n.est <- naive.est.opt(x=x, y=y,r=r, nseg = Ks[j], bdeg = bdeg,pord = pord,
                             x.grid = x)
      f.hat <- gam(y~s(x,k=(Ks[j]+bdeg), bs="ps", m=c(bdeg-1, pord)),
                   method="GCV.Cp")
      sig.hat <- sigma(f.hat)
      #n.est$sig.hat
      est.emise <-  mise.lambda.optimVB(lambda=lambda, x=x,y=y, r=r,sig=sig.hat,
                                      nseg=Ks[j], pord=pord, 
                                      bdeg=bdeg, f=f.hat$fitted.values, fr=n.est$fr.hat, p=p)
      
      keep.difmise$n[counter ] <- ns[j]
      keep.difmise$K[counter ] <-Ks[j]
      keep.difmise$exp.vardiff[counter] <- (t.emise$var-est.emise$var)^2
      keep.difmise$exp.vardiffp[counter] <- (t.emise$varp-est.emise$varp)^2
      keep.difmise$exp.biasdiff[counter] <- abs(t.emise$bias-est.emise$bias)
      keep.difmise$exp.sqdiff[counter] <- (t.emise$mise-est.emise$mise)^2
      keep.difmise$sig.hats.diff[counter] <- ((sig.hat^2) -(sig^2))^2
      keep.difmise$hh.sq[counter]<- (t.emise$HH)
      keep.difmise$hh.sqp[counter]<- (t.emise$hh)
      counter <- counter +1
      print(i)
      print("true mise")
      print(t.emise)
      print("est mise")
      print(est.emise)
      print("***************")
      print(sig.hat)
      print('^^^^^^')
      #print(t.emise)
    }
    print("*****")
    print(counter)
    print("===+++===")
    print(ns[j])
    
  }
  keep.difmise1 <- keep.difmise %>% 
    group_by(n) %>% 
    summarise(K = first(K),
              exp.vardiff = mean(exp.vardiff, na.rm = TRUE),
              exp.vardiffp = mean(exp.vardiffp, na.rm = TRUE),
              exp.biasdiff = mean(exp.biasdiff, na.rm = TRUE),
              exp.sqdiff = mean(exp.sqdiff, na.rm = TRUE),
              sig.hats.diff = mean(sig.hats.diff, na.rm = TRUE),
              hh.sq = mean(hh.sq, na.rm = TRUE),
              hh.sqp = mean(hh.sqp, na.rm = TRUE)
              
    )  
  return(list(keep.difmise1=keep.difmise1,keep.difmise=keep.difmise) )
}






## penalized spline estimate of the error variance
sigma.est <- function(x=x, y=y,r=r
                          , nseg = 35, bdeg = 4,pord = 2,
                          x.grid = x.grid){
  n<- length(x)
  initial.lambda <- 1
  searh<- optim(par = initial.lambda,fn=gcvlambda ,x=x,y=y,
                nseg=nseg, pord = pord,
                bdeg=bdeg,
                method = "L-BFGS-B",lower=0,upper=Inf)
  
  idx<- searh$par
  fr.est <- pgams(x=x,y=y, lambda = idx, r=r, x.grid = x.grid, nseg = nseg, 
                  pord = pord,
                  bdeg = bdeg)
  sig.hat <- sqrt(sum((y-fr.est$f.hat)^2)/(n-sum(diag(fr.est$A))))
  
  return(sig.hat=sig.hat)
}




sig.sim <- function(r=1, bdeg=3, pord=2,sig=0.1,nsim=5,ns){
  library(dplyr)
  #ns <- samp(range_start=100,range_end=1000,points_per_range=l)
  v <- ( (2 * (2*pord))+ 1)
  #Ks <- ceiling(15*ns^((1/v)  ))
  Ks <- floor(25* (ns^(1/(v))) )
  keep.sigmas <- data.frame(
    n = rep(NA, nsim*length(ns)),
    K = rep(NA, nsim*length(ns)),
    sigma.gcv= rep(NA, nsim*length(ns)),
    sigma.m= rep(NA, nsim*length(ns)),
    sigma.mdiff= rep(NA, nsim*length(ns)),
    sigma.gcvdiff= rep(NA, nsim*length(ns)),
    edf=rep(NA, nsim*length(ns)),
    trace =rep(NA, nsim*length(ns))
  )
  counter <- 1
  for(j in 1: length(ns)){

    set.seed(200)
    x <- seq(0, 1, length.out=ns[j])
    f <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
    if(r==1){
      fr <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
    }else if(r==2){
      fr <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
    }
    
    for (i in 1:nsim) {
      y <- f + rnorm(n=ns[j],mean=0,sd=sig)
      
      f.hat <- gam(y~s(x,k=(Ks[j]+bdeg), bs="ps", m=c(bdeg-1, pord)),
                   method="GCV.Cp")
      sig_hat <- sigma(f.hat)
      # sig_hat
      tr <- naive.est.opt(x=x, y=y,r=r
                      , nseg = Ks[j], bdeg = bdeg,pord = pord,
                      x.grid = x)
      
      
      
      keep.sigmas$n[counter ] <- ns[j]
      keep.sigmas$K[counter ] <-Ks[j]
      keep.sigmas$sigma.gcv[counter] <- sig_hat
      keep.sigmas$sigma.m[counter] <- tr$sig.hat
      keep.sigmas$sigma.mdiff[counter] <- (tr$sig.hat-sig)^2
      keep.sigmas$sigma.gcvdiff[counter] <- (sig_hat-sig)^2
      keep.sigmas$edf[counter]<- tr$edf
      keep.sigmas$trace[counter]<- tr$tr
      counter <- counter +1
      print(i)
      print("manual")
      print(tr$sig.hat)
      print("GCV")
      print(sig_hat)
    }
    print("*****")
    print(counter)
    #print(n.est$sig.hat)
  }
  keep.sigma <- keep.sigmas %>% 
    group_by(n) %>% 
    summarise(K = first(K),
              sigma.gcv = mean(sigma.gcv, na.rm = TRUE),
              sigma.m = mean(sigma.m, na.rm = TRUE),
              sigma.mdiff = mean(sigma.mdiff, na.rm = TRUE),
              sigma.gcvdiff = mean(sigma.gcvdiff, na.rm = TRUE),
              edf = mean(edf, na.rm = TRUE),
              trace = mean(trace, na.rm = TRUE)
              )  
  return(list(keep.sigma=keep.sigma,keep.sigmas=keep.sigmas) )
}



##### A more robust code
mise.diff <- function(lambda=0.1, x,y, r=1,sig=0.1,nseg=35, pord=2, 
                                bdeg=35, f=f, fr=fr,target_x=target_x)
{
  n<- length(x)
  idx <- which(x == target_x)
  # y <- f + rnorm(n=n,mean=0,sd=sig)
  f.hat <- gam(y~s(x,k=(nseg+bdeg), bs="ps", m=c(bdeg-1, pord)),
               method="GCV.Cp")
  sig.hat <- sigma(f.hat)
  naive <- naive.est.opt(x=x, y=y,r=r
                      , nseg = nseg, bdeg = bdeg,pord = pord,
                      x.grid = x)
  fr.hat <- naive$fr.hat
  BS =BS<- Bbase(x, nseg = nseg, bdeg = bdeg)
  dx <- (BS$xr- BS$xl)/nseg
  bdeg <- BS$bdeg
  knots <- BS$knots
  B <- BS$B
  v <- dim(B)[2]
  Dm <- diff(diag(v), differences = pord)
  Pm <- t(Dm) %*% Dm
  Dr <- diff(diag(v), differences =r)
  B.qr<- Bbase(x, nseg = nseg, bdeg = bdeg-r)
  Btil <-as.matrix((B.qr$B%*%Dr)/(dx^r))
  # dim(Btil)
  A.prime.lam <- Btil %*% solve((t(B)%*%B + lambda*Pm))%*%t(B) 
  # A.prime.lam.f <- as.vector(A.prime.lam%*%f)
  # dim(A.prime.lamp)
  b <- sum(A.prime.lam[idx,]^2)
  b1 <- sum(A.prime.lam^2)
  sig.sqdiff<- abs(sig^2- sig.hat^2)/n
  var.diff <- ((sig.sqdiff)*b1)
  varp.diff <- (((sig.sqdiff))*b)
  
 a1<-  sum(((A.prime.lam%*%f)^2 -(A.prime.lam%*%f.hat$fitted.values)^2)/n)
 a2<-  sum(((fr.hat)^2 -(fr)^2)/n)
 a3 <- sum(2*((A.prime.lam%*%f.hat$fitted.values)*fr.hat-(A.prime.lam%*%f)*fr )/n)
 sq.bias.diff <- abs(a1+a2+a3)
 # plot(x, (A.prime.lam%*%f)*fr)
 # plot(x, (A.prime.lam%*%f.hat$fitted.values)*fr.hat)
  return(list( mise.diff=(var.diff+sq.bias.diff), 
               var.diff=var.diff,var.diffp=varp.diff,
               sq.bias.diff=sq.bias.diff,a1=a1,a2=a2,a3=a3,
               HH = (b1), hh=(b), sig.sqdiff=sig.sqdiff ) )
}



emise.diff <- function(lambda=0.2,r=1, bdeg=3, pord=2,sig=0.1,nsim=5,
                     log.n=seq(log10(1000), log10(10000),log10(1.11) )){
  library(dplyr)
  #ns <- samp(range_start=100,range_end=1000,points_per_range=l)
  v <- ( (2 * (2*pord))+ 1)
  # Ks <- ceiling(10*ns^((1/v)  ))
  # log.n<- seq(log10(1000), log10(10000),log10(1.11) )
  ns=ceiling(10^log.n)
  rate =1/v
  Ks <-ceiling( 10^( 1.6+rate*log.n) )
  keep.difmise <- data.frame(
    n = rep(NA, nsim*length(ns)),
    K = rep(NA, nsim*length(ns)),
    var.diff= rep(NA, nsim*length(ns)),
   var.diffp= rep(NA, nsim*length(ns)),
    sq.bias.diff= rep(NA, nsim*length(ns)),
    mise.diff = rep(NA, nsim*length(ns)),
    sig.sqdiff = rep(NA, nsim*length(ns)),
    hh.sq = rep(NA, nsim*length(ns)),
    hh.sq.p = rep(NA, nsim*length(ns))
  )
  counter <- 1
  for(j in 1:length(ns)){
    #set.seed(200)
    x <- seq(0, 1, length.out=ns[j])
    f <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
    if(r==1){
      fr <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
    }else if(r==2){
      fr <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
    }
    # p <- floor(ns[j]/2)
    for (i in 1:nsim) {
      y <- f + rnorm(n=ns[j],mean=0,sd=sig)
      t.emise <- mise.diff(lambda=lambda, x,y, r=r,sig=sig,nseg=Ks[j], pord=pord, 
                         bdeg=bdeg, f=f, fr=fr)
      keep.difmise$n[counter ] <- ns[j]
      keep.difmise$K[counter ] <-Ks[j]
      keep.difmise$var.diff[counter] <- t.emise$var.diff
      keep.difmise$varp.diff[counter] <- t.emise$var.diffp
      keep.difmise$sq.bias.diff[counter] <- t.emise$sq.bias.diff
      keep.difmise$mise.diff[counter] <- t.emise$mise.diff
      keep.difmise$sig.sqdiff[counter] <- t.emise$sig.sqdiff
      keep.difmise$hh.sq[counter]<- (t.emise$HH)
      keep.difmise$hh.sqp[counter]<- (t.emise$hh)
      counter <- counter +1
      print(i)
      print("emise")
      print(t.emise)
      print("***************")
      print('^^^^^^')
      #print(t.emise)
    }
    print("*****")
    print(counter)
    print("===+++===")
    print(ns[j])
    
  }
  keep.difmise1 <- keep.difmise %>% 
    group_by(n) %>% 
    summarise(K = first(K),
              var.diff = mean(var.diff, na.rm = TRUE),
              varp.diff = mean(var.diff, na.rm = TRUE),
              sq.bias.diff = mean(sq.bias.diff, na.rm = TRUE),
              mise.diff = mean(mise.diff, na.rm = TRUE),
              sig.sqdiff = mean(sig.sqdiff, na.rm = TRUE),
              hh.sq = mean(hh.sq, na.rm = TRUE),
              hh.sqp = mean(hh.sqp, na.rm = TRUE)
              
    )  
  return(list(keep.difmise1=keep.difmise1,keep.difmise=keep.difmise) )
}





#### Study the behavior of Atilde and A#######
AtildeA <- function(x=x,y=y,lambda=0.1, r=0, x.grid, nseg=35,pord=2,bdeg=3){
  n <- length(x)
  BS<- Bbase(x, nseg = nseg, bdeg = bdeg)
  BS.grid<- Bbase(x.grid, nseg = nseg, bdeg = bdeg)
  B.grid <-BS.grid$B
  dx <- (BS$xr- BS$xl)/nseg
  knots <- BS$knots
  B <- BS$B
  ####. Cholesky Decomposition
  R <- chol(t(B)%*%B)
  ### making my D matrix
  v=dim(BS$B)[2]
  Dm <- diff(diag(v), diff = pord)
  Mv <- solve(t(R))%*%t(Dm)%*%Dm%*%solve(R)
  ## Performing svd to find the Gamma and T inverse 
  ssvd <- svd(Mv)
  Tinves <- solve(R)%*%ssvd$u
  Gamma <- diag(ssvd$d)
  #gama<-mean(diag(Gamma))
  ### Finding reparametrised parameters
  hK <- solve(diag(1, nrow = v)+lambda*Gamma)%*%(t(B%*%Tinves))
  # Thetahat <- hK%*%y
  Atrans <- t(B%*%Tinves)
  # f.hat <- A%*%y
  # 
  #theta <-1/(1+lambda*gama)*diag(1, nrow = v)%*%t(B%*%Tinves)%*%y
  ## Extract the alpha hat(original estimates)
  # alpha.hat <-  Tinves %*% Thetahat
  # fg.hat <- B.grid%*%alpha.hat
  ##### Fitting the derivative function
  if(r==0){
    Dr <- diag(v)
  }else{
    Dr <- diff(diag(v), differences = r)
    
  }
  Bqrx<- Bbase(x, nseg = nseg, bdeg = bdeg-r)
  # Bqrx.grid<- Bbase(x.grid, nseg = nseg, bdeg = bdeg-r)
  
  ## Atilde
  Atilde<- Bqrx$B%*%Dr%*%Tinves
  # Bq = (Bqrx$B%*%Dr)/(dx^r)
  # Bq.grid = (Bqrx.grid$B%*%Dr)/(dx^r)
  # K <- Bq%*%Tinves%*%hK
  # fr.hat <- Bq%*%alpha.hat
  # frg.hat<- Bq.grid%*%alpha.hat
  return(list(M=hK, Atilde=Atilde,
              AT=Atrans))
  
}

infinity_norm_matrix <- function(A) {
  max(rowSums(abs(A)))
}
ps <- function( ns, bdeg = bdeg, pord = pord, r=r, lambda = lambda){
  v <- (2 * (2 * pord)) + 1  # Smoothness parameter
  rate <- 1 / v
  c <- 1.8
  Ks <- ceiling(10^(c + rate * log10(ns)))
  print(Ks)
  ls <- sapply(ns, function(n){
    tb(n=n, bdeg = bdeg, pord = pord, r=r, lambda = lambda)
  })
  
  y1 <- unlist(ls["g",])
  y1
  y2 <- unlist(ls["gg",])
  y2
  y3 <- unlist(ls["tt",])
  y3
  y4 <- unlist(ls["ff",])
  y4
  y5 <- unlist(ls["aj",])
  y5
  y6 <- unlist(ls["aa",])
  y6
  y7 <- unlist(ls["uu",])
  y7
  y8 <- unlist(ls["utt",])
  y8
  y9 <- unlist(ls["tinv",])
  y9
  y10 <- unlist(ls["rinf",])
  y10
  y11 <- unlist(ls["rinft",])
  y11
  y12 <- unlist(ls["ll",])
  y12
  par(mfrow=c(3,4))
  plot(Ks, y1, type = "p", col="red", main = "AtildeJAT")
  plot(Ks, y2, type = "p", col="red", main ="AtildeJAT/h^r")
  plot(Ks, y3, type = "p", col="red",main ="JAT")
  plot(Ks, y4, type = "p", col="red",main ="Atilde")
  plot(Ks, y5, type = "p", col="red",main ="AtildeJ")
  plot(Ks, y6, type = "p", col="red",main ="AT")
  plot(Ks, y7, type = "p", col="red",main ="U")
  plot(Ks, y8, type = "p", col="red",main ="UT")
  plot(Ks, y9, type = "p", col="red",main ="Tinverse")
  plot(Ks, y10, type = "p", col="red",main ="Rinverse")
  plot(Ks, y11, type = "p", col="red",main ="RinverseT")
  plot(Ks, y12, type = "p", col="red",main ="TinverseT")
  # Fit linear models
  fit1 <- lm(log10(y1) ~ log10(ns))
  fit2 <- lm(log10(y2) ~ log10(ns))
  fit3 <- lm(log10(y3) ~ log10(ns))
  fit4 <- lm(log10(y4) ~ log10(ns))
  fit5 <- lm(log10(y5) ~ log10(ns))
  fit6 <- lm(log10(y6) ~ log10(ns))
  fit7 <- lm(log10(y7) ~ log10(ns))
  fit8 <- lm(log10(y8) ~ log10(ns))
  fit9 <- lm(log10(y9) ~ log10(ns))
  fit10 <- lm(log10(y10) ~ log10(ns))
  fit11 <- lm(log10(y11) ~ log10(ns))
  fit12 <- lm(log10(y12) ~ log10(ns))
  
  # Extract slopes
  slope1 <- coef(fit1)[2]
  slope2 <- coef(fit2)[2]
  slope3 <- coef(fit3)[2]
  slope4 <- coef(fit4)[2]
  slope5 <- coef(fit5)[2]
  slope6 <- coef(fit6)[2]
  slope7 <- coef(fit7)[2]
  slope8 <- coef(fit8)[2]
  slope9 <- coef(fit9)[2]
  slope10 <- coef(fit10)[2]
  slope11 <- coef(fit11)[2]
  slope12 <- coef(fit12)[2]
  # Create a data frame for the table
  slopes_table <- data.frame(
    Plot_Title = c("AtildeJAT", "AtildeJAT/h^r", "JAT", "Atilde", "AtildeJ", "AT", "U", "UT",
                   "Tinverse","Rinverse","RinverseT","TinverseT" ),
    Slope = c(slope1, slope2, slope3, slope4, slope5, slope6, slope7, slope8,
              slope9,slope10,slope11,slope12)
  )
  return(slopes_table)
}
###### Experiment
tb <-  function( n, bdeg, pord=2, r,lambda=4, c=1.4){
  v <- (2 * (2 * pord)) + 1  # Smoothness parameter
  rate <- 1 / v
  nseg <- ceiling(10^(c + rate * log10(n)))  # Number of segments
  x <- seq(0,1, length.out=n)
  BS <- Bbase(x, xl = min(x), xr = max(x), nseg = nseg, bdeg = bdeg)
  BR <- Bbase(x, xl = min(x), xr = max(x), nseg = nseg, bdeg = bdeg-r)
  B=BS$B
  dx <- BS$knots[4]-BS$knots[3] 
  sum(apply(B, 2, sum))
  apply(B, 2, sum)
  apply(B, 1, sum)
  tB <- t(B)
  sum(tB[2,])
  v=dim(BS$B)[2]
  hh <- solve(tB%*%B)
  infinity_norm_matrix(hh)
  R <- chol(t(B)%*%B)
  ccc <- solve(t(B)%*%B)
  tR <- t(R)
  infinity_norm_matrix(tR%*%R)
  infinity_norm_matrix(solve(R))
  infinity_norm_matrix(solve(tR))
  Dm <- diff(diag(v), diff = pord)
  Dr <- diff(diag(v), diff = r)
  Mv <- solve(t(R))%*%t(Dm)%*%Dm%*%solve(R)
  infinity_norm_matrix(Mv)
  
  ## Performing svd to find the Gamma and T inverse 
  ssvd <- svd(Mv)
  ee <- sort(ssvd$d, decreasing = F)
  ff<-  ee[pord+1]
  Tinves <- solve(R)%*%ssvd$u
 
    # infinity_norm_matrix(solve(R)%*%t(solve(R)) )
  rinft <- infinity_norm_matrix(t(solve(R)))
  tinv <- infinity_norm_matrix(Tinves)
  infinity_norm_matrix(ssvd$u)
  tt1 <- t(solve(R))
  ut<- t(ssvd$u)
  infinity_norm_matrix(Dr%*%Tinves)
  Gamma <- diag(ssvd$d)
  infinity_norm_matrix(Gamma)
  # lambda <- 4
  J <- solve(diag(1, nrow = v)+lambda*Gamma)
  Jt<- infinity_norm_matrix(t(J))
  Atild<- BR$B%*%Dr%*%Tinves
  A <- B%*%Tinves
  rinf <- (sum(diag(J)^2))
  ainf <- (sum(diag(J)^2))^(1/2)
    # infinity_norm_matrix(A) ## a change happened here 
  AT <- t(B%*%Tinves)
  # ff <-infinity_norm_matrix(Atild)
  aa <- infinity_norm_matrix(tB)
  tt <- infinity_norm_matrix(J%*%AT)
  gg<- infinity_norm_matrix( (Atild%*%J%*%AT)/(dx^r) )
  g<- infinity_norm_matrix( (Atild%*%J%*%AT) )
  g1<- infinity_norm_matrix( (Atild%*%AT) )
  ll<- infinity_norm_matrix(t(Tinves))
  aj <- infinity_norm_matrix((Atild%*%J))
  utt <- infinity_norm_matrix(ut)
  uu <- infinity_norm_matrix(t(ut))
  return(list(tt=tt, ainf=ainf, Jt=Jt,aa=aa, g=g, g1=g1,
              gg=gg, ff=ff, ll=ll,aj=aj, utt=utt, uu=uu, tinv=tinv,
              rinf=rinf, rinft=rinft))
}
