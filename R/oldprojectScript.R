## Load required libraries
library(parallel)
library(mgcv)
#library(tidyverse)
#library(ggplot2)

#### Making plots of the mean regression functions and their derivatives
set.seed(200)
n <- 600
x <- seq(0, 1, length.out=n)
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
fa.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
fa.pprime  <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
fb <- 4 * sin(6 * x)
fb.prime <- 24 * cos(6 * x)
fb.pprime <- -144*sin(6*x)
# Create a white background plot
par(mfrow=c(2,3))
plot(x, fa, type = "l", col = "blue", xlab = "x", ylab = paste0("f.", 'a','(x)'))
plot(x, fa.prime, type="l", col = "blue", ylab = paste0("f '.", 'a','(x)'))
plot(x, fa.pprime, type="l", col = "blue", ylab = paste0("f ''.", 'a','(x)'))

plot(x, fb, col = "red", type="l", ylab = paste0("f.", 'b','(x)'))
plot(x, fb.prime, col = "red", type="l", ylab = paste0("f '.", 'b','(x)'))
plot(x, fb.pprime, col = "red", type="l", ylab = paste0("f ''.", 'b','(x)'))


###### Functions for the estimating components #################

tpower <- function(x, t, p){
  ## Truncated p-th power function
  ## x is the values to be evaluated at knot t
  ## p is the degree of the plynomial
  
  
  ((x-t)^p)*(x>t)
}

## Creating a B-Spline basis of degree "deg"
bbase <- function(x, xl, xr, ndx, deg){
  # Construct a B-spline basis of degree ’deg’
  ## xl is the smallest value of x
  ## xr is the largest value of x
  ## ndx is length of x
  dx <- (xr - xl) / ndx
  knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
  P <- outer(x, knots, tpower, deg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
  B <- (-1)^(deg + 1) * P %*% t(D)
  return(list(B=B, knots=knots, dx=dx, deg=deg))
}

#### B-spline for new x(x.grid)
bbase.grid <- function(x,dx,knots,deg){
  P <- outer(x,knots,tpower,deg)
  n <- dim(P)[2]
  D <- diff(diag(n),diff=deg+1)/(gamma(deg+1)*dx^deg)
  B <- (-1)^(deg+1)*P %*% t(D)
  return(list(B=B,dx=dx,knots=knots,deg=deg))
}


### Function to estimate the f and its derivatives

Pgam <- function(x, y, B, r, lambda, x.grid){  
  dx <- B$dx
  deg <- B$deg
  pen.order <- deg-1
  knots <- B$knots  
  B <- B$B
  v <- dim(B)[2]
  n <- dim(B)[1]
  Dm <- diff(diag(v), diff = pen.order)
  Pm <- t(Dm) %*% Dm
  alph.hat <- solve(t(B)%*%B + lambda*Pm, t(B)%*%y)
  B.grid <- bbase.grid(x= x.grid, dx=dx, knots = knots, deg)
  ## Computing the estimated function
  f.hat <- B.grid$B %*% alph.hat
  ## computing the r-th derivative of the f
  Dr <- diff(diag(v), diff=r)
  B.qr <- bbase.grid(x= x.grid, dx =dx, knots =knots[(r+1):(length(knots)-r)], deg= deg-r)
  fr.hat <- B.qr$B%*%Dr%*%alph.hat/(dx^r)
  return(data.frame(x.grid,f.hat, fr.hat)) 
}


##Optimization

## Function to optimize
pluginmise <- function(lambda, x, B, r,sig, f.grid,fr.grid){
  dx <- B$dx
  deg <- B$deg
  pen.order <- deg-1
  knots <- B$knots
  B <- B$B
  v <- dim(B)[2]
  Dm <- diff(diag(v), diff = pen.order)
  Pm <- t(Dm) %*% Dm
  Dr <- diff(diag(v), diff=r)
  B.qr <- bbase.grid(x= x, dx =dx, knots =knots[(r+1):(length(knots)-r)], deg= deg-r)
  
  Btil <-as.matrix(B.qr$B %*%Dr/(dx^r))
  
  A.prime.lam <- Btil %*% solve(t(B)%*%B + lambda*Pm)%*%t(B) 
  
  variance <- (sig^2)*sum(A.prime.lam^2)/length(x.grid)
  sq.bias <- sum(as.vector(A.prime.lam%*%f.grid - fr.grid)^2)/length(x.grid)
  
  MISE.approx <- variance+sq.bias
  return(MISE.approx)
}



fit.f.naive.fr <- function(x,y,deg,r,k,x.grid){
  ### This function fits the naive estimator on a grid
  order <- deg-1
  ndx <- k-deg
  fit <- gam(y~s(x,k=k,bs="ps",m=c(order,order-1)),method="GCV.Cp")
  beta.hat <- fit$coef # get cofficients
  # get design matrix from gam fit
  mat <- predict.gam(fit, type = "lpmatrix")
  # get knots and diff between them
  knots <- fit$smooth[[1]]$knots
  dx <- diff(knots)[1] 
  # make Bpline Basis
  Bx <- (bbase.grid(x,dx,knots,deg))$B
  alpha.hat <- as.vector(solve(t(Bx)%*%Bx)%*%t(Bx)%*%mat%*%beta.hat)
  D=as.matrix(fit$smooth[[1]]$D)
  p = t(D)%*%D
  lambdaGCV1= (t(Bx)%*%y- t(Bx)%*%Bx%*%alpha.hat)/(p%*% alpha.hat)
  lambdagams1 =lambdaGCV1[1]
  # alpha <- fit$smooth[[1]]$S.scale
  # lambdaGCV2 <- fit$sp /alpha
  # lambdagams2 =lambdaGCV2[1]
  # 
  ##  estimate of f on grid
  B.grid <- bbase.grid(x.grid, dx,knots,deg)
  f.hat <- as.vector(B.grid$B %*% alpha.hat)
  v <- dim(B.grid$B)[2]
  ## naive r-th derivative estimate on grid
  Dr<- diff(diag(v), diff=r)
  B.qr <- bbase.grid(x= x.grid, dx =dx, knots =knots[(r+1):(length(knots)-r)], deg= deg-r)
  fr.hat <- B.qr$B%*%Dr%*%alpha.hat/(dx^r)
  return(data.frame(f.hat=f.hat, naive.fr.hat =fr.hat, est=lambdagams1))
}



est.f.fr.plug.in <- function(x, y, deg, r, k,x.grid){
  
  xr <- max(x)-0.01
  xl<- min(x)+0.01
  order <- deg-1
  ndx <- k-deg
  f.hat <- gam(y~s(x,k=k,bs="ps",m=c(order,order-1)),method="GCV.Cp")
  ## Using estimated Sigma.
  sig.hat <- sigma(f.hat)
  beta.hat <- f.hat$coef # get cofficients
  # get design matrix from gam fit
  mat <- predict.gam(f.hat, type = "lpmatrix")
  # get knots and diff between them
  knots <- f.hat$smooth[[1]]$knots
  dx <- diff(knots)[1] 
  # make Bpline Basis
  Bx <- (bbase.grid(x,dx,knots,deg))
  # get alpha.hat for Bx (change of basis)
  alpha.hat <- as.vector(solve(t(Bx$B)%*%Bx$B)%*%t(Bx$B)%*%mat%*%beta.hat)
  ##########Fitting the naive estimator using alpha hat##########

  v <- dim(Bx$B)[2]
  ## f estimator using GCV selected lambda
  B.grid <- bbase.grid(x= x.grid, dx=dx, knots = knots, deg)
  f.hat <- as.vector(B.grid$B %*% alpha.hat)
  ## est using the data
  Dr <- diff(diag(v), diff=r)
  B.qr <- bbase.grid(x= x.grid, dx =dx, knots =knots[(r+1):(length(knots)-r)], deg= deg-r)  
  fr.hatnaive <- as.vector(B.qr$B%*%Dr%*%alpha.hat/(dx^r))
  
  initial.lambda <- 10
  estimate <- optim(par = initial.lambda,fn=pluginmise, 
                    x=x.grid, B=B.grid, r=r ,sig = sig.hat, f.grid=f.hat, fr.grid=fr.hatnaive,
                    method = "L-BFGS-B",lower=0,upper=Inf)
  
  
  lamdaN <- (estimate$par)
  #BS <- bbase(x=x, xl=xl,xr=xr, ndx=ndx, deg = deg)
  modelN <- Pgam(x=x,y=y,B=Bx,r=r,lambda = lamdaN, x.grid = x.grid)
  
  
  fr.hat.plug.in <- modelN$fr.hat
  plugin.f <- modelN$f.hat
  
  return(data.frame(x.grid=x.grid,plugin.f=plugin.f, fr.hat.plug.in =fr.hat.plug.in,est.lambda= lamdaN, sig.hat=sig.hat)) 
}



####Function to estimate trueOpt

trueopt.est <- function(x=x,y,B, x.grid,B.grid, r,sig,f.grid,fp.grid){
  ##Using the true functions
  initial_lambda <- 0.03
  Optimal_lambdaT <- optim(par = initial_lambda, fn= pluginmise, x=x.grid, B=B.grid, r=r,sig = sig, f.grid=f.grid, fr.grid=fp.grid,
                           method = "L-BFGS-B", lower = 0, upper = Inf)
  
  ## Choosing small lambda to estimate the r derivative
  lamda <- Optimal_lambdaT$par
  modelT <- Pgam(x=x,y=y,B=B,r=r,lambda = lamda, x.grid = x.grid)
  f.true <- modelT$f.hat
  fr.true <- modelT$fr.hat
  return(data.frame(x.grid=x.grid,f.true=f.true, fr.true=fr.true, est.lambda=lamda))
  
}


#### Function to estimate the Iterated

Iterated.fit.est <-function(N, x,y, x.grid, B, B.grid, r, sig, plugin.f, plugin.fr) {
  
  for(k in 1:N){
    #plugin.fr.1 <-plugin.fr
    #BS1=bbase(x=x.grid, xl=xl,xr=xr, ndx=ndx, deg = deg)
    initial_lambda <- 0.03
    Optimal_lambdaI <- optim(par = initial_lambda,fn= pluginmise, x=x.grid, B=B.grid, r=r,sig = sig, f.grid=plugin.f, fr.grid=plugin.fr,
                             method = "L-BFGS-B", lower = 0, upper = Inf)
    
    ## Choosing small lambda to estimate the r derivative
    lamdaI <- Optimal_lambdaI$par
    #BS=bbase(x=x, xl=xl,xr=xr, ndx=ndx, deg = deg)
    modelI <- Pgam(x=x,y=y,B=B,r=r,lambda = lamdaI, x.grid = x.grid)
    
    
    improv <- mean((plugin.fr-modelI$fr.hat )^2)
    print("******")
    print(k)
    print(improv)
    plugin.f <- modelI$f.hat
    plugin.fr<- modelI$fr.hat
    if(k ==N-1){
      plugin.frN1 <- plugin.fr
    }
  }
  
  
  Iterated.f <- plugin.f
  Iterated.fr <- plugin.fr
  
  
  return(data.frame(x.grid=x.grid, Iterated.f=Iterated.f,plugin.frN1=plugin.frN1,Iterated.fr=Iterated.fr, est.lamda=lamdaI ))
  
}






## Data
set.seed(200)
n <- 600
sig= 1
x <- seq(from=0,to=1,length=n)
x.grid <- seq(from=min(x),to=max(x),length=100)

y <- 32*exp(-8*(1-2*x)^2)*(1-2*x)+sig*rnorm(n)
f <-32*exp(-8*(1-2*x)^2)*(1-2*x)
fp <- (4096*x^2-4096*x+960)*exp(-8*(1-2*x)^2)
f.grid <-32*exp(-8*(1-2*x.grid)^2)*(1-2*x.grid)
fp.grid <- (4096*x.grid^2-4096*x.grid+960)*exp(-8*(1-2*x.grid)^2)
r <- 1
deg <- 3
k <- 35
ndx <- k-deg
xr <- max(x)
xl<- min(x)
# test the functions

BS <- bbase(x, xl=xl,xr=xr, ndx=ndx, deg=deg)

pointwiseCI <- function(deriv.estimate,deriv.order, sig, lambda, x.grid,BS){
  B<- BS$B
  deg<- BS$deg
  knots<- BS$knots
  dx<- BS$dx
  pen.order <- deg-1
  v <- dim(B)[2]
  r<-deriv.order
  Dm <- diff(diag(v), diff = pen.order)
  Pm <- t(Dm) %*% Dm
  dim(B)
  inv.part<- solve( t(B)%*%B+lambda*Pm)
  B.grid <- bbase.grid(x.grid,dx,knots,deg)
  Dr <- diff(diag(v), diff=r)
  B.qr <- bbase.grid(x= x.grid, dx =dx, knots =knots[(r+1):(length(knots)-r)], deg= deg-r)
  est.var <- (sig^2)*(B.qr$B%*%Dr)%*%inv.part %*%t((B.grid$B)) %*%t((B.qr$B%*%Dr)%*%inv.part %*%t((B.grid$B))) /(dx^(2*r) )
  est.st.dev<- sqrt(diag(est.var))
  
  pointwiseCI.lower <- deriv.estimate-2*est.st.dev
  pointwiseCI.upper <- deriv.estimate+2*est.st.dev
  
  return(list(Estimate= deriv.estimate, CI.lower=pointwiseCI.lower,CI.upper=pointwiseCI.upper))
}

B.grid <- bbase.grid(x.grid,BS$dx,BS$knots,deg)
dim(BS$B)

## Pgam function
est <- Pgam(x, y, B=BS, r=1, lambda=0.9, x.grid)
plot(est$x.grid, est$f.hat, type="l")
plot(est$x.grid, est$fr.hat, type="l")



cnf <- pointwiseCI(deriv.estimate =est$fr.hat, deriv.order = 1, 
            sig = 1,lambda=0.9, x.grid=x.grid, BS=BS)

plot(x.grid, cnf$Estimate, type = "l")
lines(x.grid, cnf$CI.lower, col="red", lty=2)
lines(x.grid, cnf$CI.upper, col="red",lty=2)
## MISE function
a <- pluginmise(lambda=0.2, x=x.grid, r=r,B=B.grid, sig = sig, f.grid=f.grid, fr.grid=fp.grid )

l.grid <- seq(from=exp(-1),to=exp(5),length=101)
keep <- rep(NA,101)

for (i in 1:101){
  keep[i] <- pluginmise(lambda = l.grid[i], x=x.grid, r=r,B=B.grid, sig = sig, f.grid=f.grid, fr.grid=fp.grid )
}
plot(l.grid,keep,type="l", col="red")


naive.fits <- fit.f.naive.fr(x=x,y=y,deg,r,x.grid,k=k )

###Plotting f estimate
plot(x.grid,naive.fits$f.hat,type="l", col="red")
lines(x,f,pch=16,cex=.5)
plot(x.grid,naive.fits$naive.fr.hat,type="l")
lines(x.grid,fp.grid,pch=16,cex=.5, col="red")


## Testing the plugin estimator
fr.hat.plug.in <- est.f.fr.plug.in(x, y, deg=3, r=1, k=k,x.grid=x.grid)
plot(x.grid,fp.grid, type = "l" )
lines(fr.hat.plug.in$x.grid,fr.hat.plug.in$fr.hat.plug.in,col="blue")


## Testing the oracle estimate function
true.est <- trueopt.est(x=x,y=y,B=BS, x.grid=x.grid,B.grid=B.grid, r=r,sig=sig,f.grid=f.grid,fp.grid=fp.grid)
plot(x.grid,fp.grid, type = "l" )
lines(true.est$x.grid,true.est$fr.true,col="blue")


Iterated.est <- Iterated.fit.est(N=100, x=x,y=y, x.grid=x.grid, B=BS,
                                 B.grid=B.grid, r=r, sig=fr.hat.plug.in$sig.hat[1],
                                 plugin.f= fr.hat.plug.in$plugin.f, plugin.fr =fr.hat.plug.in$fr.hat.plug.in)

plot(x.grid,fp.grid, type = "l" )
lines(Iterated.est$x.grid,Iterated.est$Iterated.fr,col="blue")


###### Function for simulations
est.fun <- function(n, MCN, n.grid){   
  True.est <- list()
  plugin.est <- list()
  naive.est <- list()
  Iterated.est <- list()
  x <- seq(from=0,to=1,length=n)
  x.grid <- seq(from=min(x),to=max(x),length=n.grid)
  f<-32*exp(-8*(1-2*x)^2)*(1-2*x)
  f.grid <-32*exp(-8*(1-2*x.grid)^2)*(1-2*x.grid)
  fp.grid <- (4096*x.grid^2-4096*x.grid+960)*exp(-8*(1-2*x.grid)^2)
  rg= range(f)[2]-range(f)[1]
  fp <- (4096*x^2-4096*x+960)*exp(-8*(1-2*x)^2)
  
  
  for (i in 1:MCN) {
    
    print(i)
    sig <- 0.3
    y <- f+sig*rnorm(n)
    ## Setting params
    xl<-min(x)
    xr<-max(x)
    r <- 1
    k <- 25
    deg <- 3
    ndx <- k-deg
    
    BS <- bbase(x, xl=xl, xr= xr, ndx=ndx, deg=deg)
    B.grid <- bbase.grid(x.grid,BS$dx,BS$knots,deg)
    ########## Naive estimation- Derivative and mean function
    # fit Bspline
    Naive.est =fit.f.naive.fr(x=x,y=y,deg=deg,r=r,x.grid,k=k) 
    Naive <- Naive.est$naive.fr.hat
    
    plot(x.grid, fp.grid, type = "l")
     lines(x.grid,  Naive.est$naive.fr.hat )
        
    naive.est[[i]]<-Naive
    ### Estimating the TrueOpt
    True.Opt <- trueopt.est(x=x,y=y, B=BS, x.grid=x.grid, B.grid =B.grid, r=r, sig=sig, f.grid= f.grid, fp.grid=fp.grid)
    
    plot(x.grid, fp.grid, type = "l")
    lines(x.grid, True.Opt$fr.true )
    True.est[[i]] <-  True.Opt$fr.true ###### MC results for fr(opt, true f and fr)
    
    
    ###. Estimating the Plugin 
    plugin <-  est.f.fr.plug.in(x, y, deg=deg, r=r, k=k,x.grid=x.grid)
    
    plugin.1 <- plugin$fr.hat.plug.in
    # 
    plot(x.grid, fp.grid, type = "l")
    lines(x.grid, plugin$fr.hat.plug.in )
    plugin.est[[i]]<- plugin$fr.hat.plug.in
    
    ### Iterated Estimator 
    #### this is obtained by iterating the plug-in estimator hopefully to get
    ## a better estimator than the plug-in.
    Iterated <- Iterated.fit.est(N=10, x=x,y=y, x.grid=x.grid, B=BS,
                                     B.grid=B.grid, r=r, sig=plugin$sig.hat[1],
                                     plugin.f= plugin$plugin.f, plugin.fr =plugin.1)
    
    
    Iterated.est[[i]] <- Iterated$Iterated.fr 
    
    
    plot(x.grid, fp.grid, type = "l")
    lines(x.grid, Iterated.est$Iterated.fr, col="orange")
    lines(x.grid, Iterated.est$plugin.frN1, col="blue")
  }
  
  
  return( list(Naive =naive.est, plugin=plugin.est, True=True.est,Iterated =Iterated.est, x.grid=x.grid, fp.grid=fp.grid)) 
  
}

  est.fun(n=50, MCN=10, n.grid = 100)

mise.est <- function(t, est.fun=est.fun,n, MCN,n.grid){
  store.true = c()
  store.naive = c()
  store.plugin = c()
  store.iterated = c()
  for(i in 1:t){
    res.est<- est.fun(n=n, MCN=MCN, n.grid = n.grid)
    fp.grid<- res.est$fp.grid
     lss.n<- res.est[["Naive"]]
     dfa.n <- do.call(cbind, lss.n)
     colnames(dfa.n) <- paste0("Run_", 1:MCN) 
     mise.naive <- mean(apply(dfa.n, 2, function(x) mean((x-fp.grid)^2)))
     store.naive[i] <-mise.naive 
     
     lss.p<- res.est$plugin 
     dfa.p <- do.call(cbind, lss.p)
     colnames(dfa.p) <- paste0("Run_", 1:MCN) 
     mise.plugin <- mean(apply(dfa.p, 2, function(x) mean((x-fp.grid)^2)))
     store.plugin[i] <-mise.plugin 
     
     lss.t<- res.est[["True"]]
     dfa.t <- do.call(cbind, lss.t)
     colnames(dfa.t) <- paste0("Run_", 1:MCN) 
     mise.true <- mean(apply(dfa.t, 2, function(x) mean((x-fp.grid)^2)))
     store.true[i] <-mise.true 
     
     lss.I<- res.est[["Iterated"]]
     dfa.I <- do.call(cbind, lss.I)
     colnames(dfa.I) <- paste0("Run_", 1:MCN) 
     mise.iterated <- mean(apply(dfa.I, 2, function(x) mean((x-fp.grid)^2)))
     store.iterated[i] <-mise.iterated 
  }
  mise.dat<- data.frame(oracle=store.true, naive=store.naive, plugin=store.plugin, iterated= store.iterated   )
   return(mise.dat)  
}

 mcn.est<-  mise.est(t=30, est.fun=est.fun,n=100, MCN=10,n.grid=50)
 mcn.est
 save( mcn.est, file = "mise.RData")
 mcn.est<- get(load("mise.RData"))
 mcn.mise <- mcn.est %>% pivot_longer(cols = everything(), names_to = "method", values_to = "MISE")
 mcn.mise %>% ggplot(aes(y=log(MISE), fill= method))+geom_boxplot()+
   scale_x_discrete(breaks = NULL)

 boxplot(log(mcn.mise$MISE)~mcn.mise$method, col=c(5,6,10,4), axes=T, ylab = "log(mise)", xlab = "method")



 ############################ AN EXAMPLE ##########################################
 ## Implementing the plugin estimator on real data
 ## An example using Raman shift.
 
 ### example data
 ### Create the R object  BCE1  by reading in processed Raman spectrum data (unoriented 532 nm)
 ### from   http://rruff.info/bastnasite/display=default/R050409
 ###			(Delete first and last rows with  0  in second column)
 ###			This is from the first sample of cerium bastnasite
 # Load necessary library
 library(readr)
 
 # Read the data from the new URL
 data_url <- "https://rruff.info/repository/sample_child_record_raman/by_minerals/BastnasiteCe__R050409__Raman__532__0__unoriented__Raman_Data_Processed__23513.txt"
 raw_data <- read.table(data_url, header = F, sep = ",")
 colnames(raw_data)<- c("Raman.Shift", "Intensity")
 head(raw_data)
 
 # Delete first and last rows with 0 in the second column
 raw_data_clean <- raw_data[raw_data$Intensity != 0, ]
 head(raw_data_clean )
 tail(raw_data_clean )
 # Create the R object BCE1
 BCE1 <- raw_data_clean
 
 # View the first few rows of BCE1
 head(BCE1)
 par(mfrow=c(3,1))
 plot(BCE1$Raman.Shift, BCE1$Intensity )
 
 ### Create the R object  BCE5  by reading in processed Raman spectrum data (unoriented 532 nm)
 ### from   http://rruff.info/bastnasite/display=default/R060737
 ###			(Delete first and last rows with  0  in second column)
 ###			This is from the second sample of cerium bastnasite
 
 data.url <- "https://rruff.info/repository/sample_child_record_raman/by_minerals/BastnasiteCe__R060737__Raman__532__0__unoriented__Raman_Data_Processed__12623.txt"
 
 raw_data.1 <- read.table(data.url, header = F, sep = ",")
 colnames(raw_data.1)<- c("Raman.Shift", "Intensity")
 head(raw_data.1)
 
 # Delete first and last rows with 0 in the second column
 raw_data_clean.1 <- raw_data.1[raw_data.1$Intensity != 0, ]
 head(raw_data_clean.1 )
 tail(raw_data_clean.1 )
 # Create the R object BCE1
 BCE5 <- raw_data_clean.1
 
 # View the first few rows of BCE1
 head(BCE5)
 
 plot(BCE5$Raman.Shift, BCE5$Intensity )
 
 
 
 ### Create the R object  BLA   by reading in processed Raman spectrum data (unoriented 532 nm)
 ### from   http://rruff.info/bastnasite/display=default/R070142
 ###			(Delete first and last rows with  0  in second column)
 ###			This is from the sample of lanthanum bastnasite
 
 
 data.url.1 <- "https://rruff.info/repository/sample_child_record_raman/by_minerals/BastnasiteLa__R070142__Raman__532__0__unoriented__Raman_Data_Processed__19566.txt"
 
 raw.data.1 <- read.table(data.url.1, header = F, sep = ",")
 colnames(raw.data.1)<- c("Raman.Shift", "Intensity")
 head(raw.data.1)
 
 # Delete first and last rows with 0 in the second column
 raw_data.clean.1 <- raw.data.1[raw.data.1$Intensity != 0, ]
 head(raw_data.clean.1 )
 tail(raw_data.clean.1 )
 # Create the R object BCE1
 BLA <- raw_data.clean.1
 
 # View the first few rows of BCE1
 head(BLA)
 
 plot(BLA$Raman.Shift, BLA$Intensity )
 
 
 
 BCE1r <- BCE1[ ( BCE1[,1]> 800 ) & (BCE1[,1] < 1000) , ]
 #plot(BCE1r$Raman.Shift, BCE1r$Intensity)
 BCE5r <- BCE5[ ( BCE5[,1]> 800 ) & (BCE5[,1] < 1000) , ]
 BLAr <-  BLA[ ( BLA[,1]> 800 ) & (BLA[,1] < 1000) , ]

 # 
 # 
 # xBCE1r <-  ( BCE1r[,1] - 900 )  / 100
 # xBCE5r <-  ( BCE5r[,1] - 900 )  / 100
 # xBLAr  <-  ( BLAr[,1]  - 900 )  / 100
 # 
 n <- 415
 xfit <- 2*((1:n)-.5)/n-1
 xfit.1 <- ((0:n))/(n)
 
 
 xBCE1r1<- (BCE1r[,1]- min(BCE1r[,1]) )/ (max(BCE1r[,1])-min(BCE1r[,1]))
 xBCE5r1<- (BCE5r[,1]- min(BCE5r[,1]) )/ (max(BCE5r[,1])-min(BCE5r[,1]))
 xBLAr1<- (BLAr[,1]- min(BLAr[,1]) )/ (max(BLAr[,1])-min(BLAr[,1]))

 # par(mfrow=c(3,1))
 plot(xBCE1r1, BCE1r$Intensity )
 
 
 plot(xBCE5r1, BCE5r$Intensity )
 plot(xBLAr1, BLAr$Intensity )

 y<- matrix(0,n,3)
 y[,1] <-  BCE1r[,2]
 y[,3] <-  BCE5r[,2]
 y[,2] <-  BLAr[,2]
 
 #RamanData <- data.frame()
 
 # plot(1,type="n",ylim=c(0,800),xlim=c(800,1035),xlab="Raman Shift",ylab="Intensity",cex.axis=0.9)
 # matplot(xfit.1, y, type = "l", lwd = c(2,2,2), col = c(2,3,4))
 # legend("topright", c("CeCO_3_F","LaCO_3_F","to be inferred"),col=c(2,3,4),pch=c(2,3,4))
 # title("(a) Experimental data")
 
 
 
 
 #### Penalized derivative estimates for the functions.
 
 
 plot( BCE1r$Raman.Shift,  BCE1r$Intensity)
 plot( BCE5r$Raman.Shift,  BCE5r$Intensity)
 plot( BLAr$Raman.Shift,  BLAr$Intensity)
 
 deg<-3
 k=35
 ndx=k-deg

 x.grid <- xfit.1
 ## Estimate the first derivative for the first sample of CeCO3F
 est.deriv.1BCE1 <- est.f.fr.plug.in(x=xBCE1r1, y=y[,1], deg = deg, r=1, k=k, x.grid =x.grid)
 naive.BCE1.fit <- fit.f.naive.fr(x=xBCE1r1, y=y[,1], deg=deg, r=1, k=k, x.grid = x.grid)

 ## plot of the function itself
plot(x.grid,  naive.BCE1.fit$f.hat)
points(xBCE1r1, y[,1])

## plot of the derivative
plot(x.grid,  naive.BCE1.fit$naive.fr.hat, type = "l", col="blue")
 
 est.deriv.1BCE1$sig.hat
 est.deriv.1BCE1$est.lambda
## using the Iterated


B <- bbase(x=xBCE1r1 ,xl=min(xBCE1r1), xr=max(xBCE1r1), ndx = ndx, deg = deg)
#B.grid <- bbase(x=xfit.1 ,xl=min(xfit.1), xr=max(xfit.1), ndx = 6, deg = 3)
B.grid <- bbase.grid(x.grid,B$dx,B$knots,deg)
#Iterated.fit.est()

tail(B.grid$B)
head(B.grid$B)
matplot(xfit.1, B.grid$B, type = "l")
Ite.BCE1.est  <- Iterated.fit.est(N=25, x=xBCE1r1, y=y[,1],B=B,r=1, B.grid = B.grid,
                                 plugin.f =est.deriv.1BCE1$plugin.f,
                                 plugin.fr =est.deriv.1BCE1$fr.hat.plug.in,
                                 x.grid =est.deriv.1BCE1$x.grid,sig = est.deriv.1BCE1$sig.hat[1] )
#Ite.BCE1.est$

CI.BCE1<- pointwiseCI(deriv.estimate =Ite.BCE1.est$Iterated.fr, deriv.order = 1, 
            sig =est.deriv.1BCE1$sig.hat[1] ,lambda = Ite.BCE1.est$est.lamda[1], 
            x.grid=x.grid, BS=B)




# plot(xBCE1r1, y[,1])










## Estimate the first derivative for the first sample of LaCO3F
est.deriv.1BLA <- est.f.fr.plug.in(x=xBLAr1, y=y[,2], deg = deg, r=1, k=k, x.grid = x.grid)
naive.BLA.fit <- fit.f.naive.fr(x=xBLAr1, y=y[,2], deg = deg, r=1, k=k, x.grid = x.grid)



est.deriv.1BLA$sig.hat

B <- bbase(x=xBLAr1 ,xl=min(xBLAr1), xr=max(xBLAr1), ndx = ndx, deg = deg)
B.grid <- bbase.grid(x.grid,B$dx,B$knots,deg)
#Iterated.fit.est()

Ite.BLA.est <- Iterated.fit.est(N=50, x=xBLAr1, y=y[,2],B=B,r=1, B.grid = B.grid,
                                 plugin.f =est.deriv.1BLA$plugin.f,
                                 plugin.fr =est.deriv.1BLA$fr.hat.plug.in,
                                 x.grid =est.deriv.1BLA$x.grid,sig = est.deriv.1BLA$sig.hat[1] )


CI.BLA<- pointwiseCI(deriv.estimate =Ite.BLA.est$Iterated.fr, deriv.order = 1, 
                      sig =est.deriv.1BLA$sig.hat[1] ,lambda = Ite.BLA.est$est.lamda[1], 
                      x.grid=x.grid, BS=B)










# plot(xBLAr1, y[,2])
# lines(x.grid,Ite.BLA.est$Iterated.f, lwd=3, col="red")




################# Estimate the first derivative for the to be inferred sample###################

est.deriv.1BCE5 <- est.f.fr.plug.in(x=xBCE5r1, y=y[,3], deg = deg, r=1, k=k, x.grid = x.grid)
naive.BCE5.fit <- fit.f.naive.fr(x=xBCE5r1, y=y[,3], deg = deg, r=1, k=k, x.grid = x.grid)

est.deriv.1BCE5$sig.hat

B <- bbase(x=xBCE5r1 ,xl=min(xBCE5r1), xr=max(xBCE5r1), ndx = ndx, deg = deg)
B.grid <- bbase.grid(x.grid,B$dx,B$knots,deg)
#Iterated.fit.est()

Ite.BCE5.est <- Iterated.fit.est(N=50, x=xBCE5r1, y=y[,3],B=B,r=1, B.grid = B.grid,
                                plugin.f =est.deriv.1BCE5$plugin.f,
                                plugin.fr =est.deriv.1BCE5$fr.hat.plug.in,
                                x.grid =est.deriv.1BCE5$x.grid,
                                sig = est.deriv.1BCE5$sig.hat[1] )


CI.BCE5<- pointwiseCI(deriv.estimate =Ite.BCE5.est$Iterated.fr, deriv.order = 1, 
                     sig =est.deriv.1BCE5$sig.hat[1] ,lambda = Ite.BCE5.est$est.lamda[1], 
                     x.grid=x.grid, BS=B)
########################################################################################





########### ploting the mean functions ################################################

par(mfrow=c(1,3))
## plot of original function comparing naive and iterated
plot(x.grid,Ite.BCE1.est$Iterated.f, lwd=3, col="red", ylim = c(0,400), type = "l",ylab = "BCE1 Intensity")
lines(x.grid,  naive.BCE1.fit$f.hat, lwd=3, col="blue")
lines(x.grid, est.deriv.1BCE1$plugin.f, lwd=3, col="green")
points(xBCE1r1, y[,1])
legend("topleft", legend = c("iterated", " naive","plugin", "data"), col=c("red", "blue","green", "black"),
       lty = 2)


## plot of original function comparing naive ,plugin and iterated
plot(x.grid,Ite.BLA.est$Iterated.f, lwd=3, col="red", ylim = c(0,400), type = "l",ylab = "BLA Intensity")
lines(x.grid,  naive.BLA.fit$f.hat, lwd=3, col="blue")
lines(x.grid, est.deriv.1BLA$plugin.f, lwd=3, col="green")
points(xBLAr1, y[,2])
legend("topleft", legend = c("iterated", " naive","plugin", "data"), col=c("red", "blue","green", "black"),
       lty = 2)


## plot of original function comparing naive ,plugin and iterated
plot(x.grid,Ite.BCE5.est$Iterated.f, lwd=3, col="red", ylim = c(0,700), type = "l",ylab = "BCE5 Intensity")
lines(x.grid,  naive.BCE5.fit$f.hat, lwd=3, col="blue")
lines(x.grid, est.deriv.1BCE5$plugin.f, lwd=3, col="green")
points(xBCE5r1, y[,3])
 legend("topleft", legend = c("iterated", " naive","plugin", "data"), col=c("red", "blue","green", "black"),
       lty = 2)
#########################################################################################################


 
 
 ########## Derivatives plots ##########################################################################
 
 par(mfrow=c(1,3))
 ## plot of derivatives function comparing naive and iterated
 plot(x.grid,Ite.BCE1.est$Iterated.fr, lwd=2, col="red", ylim = c(-5000,4000), type = "l", ylab = "BCE1 derivative Intensity")
 lines(x.grid,  naive.BCE1.fit$naive.fr.hat, lwd=2, col="blue")
 lines(x.grid, est.deriv.1BCE1$fr.hat.plug.in, lwd=2, col="green")
 #points(xBCE1r1, y[,1])
 legend("topleft", legend = c("iterated", " naive","plugin"), col=c("red", "blue","green"),
        lty = 2)
 
 
 ## plot of original function comparing naive and iterated
 plot(x.grid,Ite.BLA.est$Iterated.fr, lwd=2, col="red", ylim = c(-5000,4000), type = "l", ylab = "BLA derivative  Intensity")
 lines(x.grid,  naive.BLA.fit$naive.fr.hat, lwd=2, col="blue")
 lines(x.grid, est.deriv.1BLA$fr.hat.plug.in, lwd=2, col="green")
 #points(xBCE1r1, y[,1])
 legend("topleft", legend = c("iterated", " naive","plugin"), col=c("red", "blue","green"),
        lty = 2)
 

## plot of original function comparing naive and iterated
plot(x.grid,Ite.BCE5.est$Iterated.fr, lwd=2, col="red", ylim = c(-5000,5000), type = "l", ylab = "BCE5 derivative Intensity")
lines(x.grid,  naive.BCE5.fit$naive.fr.hat, lwd=2, col="blue")
lines(x.grid, est.deriv.1BCE5$fr.hat.plug.in, lwd=2, col="green")
#points(xBCE1r1, y[,1])
legend("topleft", legend = c("iterated", " naive","plugin"), col=c("red", "blue","green"),
       lty = 2)

##################################################################################################



par(mfrow=c(3,1))
plot(x.grid,CI.BCE1$Estimate, type = "l")
lines(x.grid,CI.BCE1$CI.lower, type = "l", col="red",lty=2)
lines(x.grid,CI.BCE1$CI.upper, type = "l", col="red",lty=2)

plot(x.grid,CI.BLA$Estimate, type = "l")
lines(x.grid,CI.BLA$CI.lower, type = "l", col="red",lty=2)
lines(x.grid,CI.BLA$CI.upper, type = "l", col="red",lty=2)

plot(x.grid,CI.BCE5$Estimate, type = "l")
lines(x.grid,CI.BCE5$CI.lower, type = "l", col="red",lty=2)
lines(x.grid,CI.BCE5$CI.upper, type = "l", col="red",lty=2)


# par(mfrow=c(1,2))
# plot(x.grid,est.deriv.1BCE5$fr.hat.plug.in, type = "l")
# lines(Ite.BCE5.est$x.grid,Ite.BCE5.est$Iterated.fr, type = "l", col="red")
# 
# 
# plot(xBCE5r1, y[,3])
# lines(x.grid,Ite.BCE5.est$Iterated.f, lwd=3, col="red")
# 




y.hat<- matrix(0,length(x.grid),3)
y.hat[,1]<-Ite.BCE1.est$Iterated.fr
y.hat[,2]<- Ite.BLA.est$Iterated.fr
y.hat[,3]<-Ite.BCE5.est$Iterated.fr


par(mfrow=c(1,2))

plot(min(BCE1r[,1])+xBCE1r1*(max(BCE1r[,1])-min(BCE1r[,1])), y[,1],type="n",ylim=c(0,800),xlim=c(800,1035),xlab="Raman Shift",ylab="Intensity",cex.axis=0.9)
points(min(BCE1r[,1])+xBCE1r1*(max(BCE1r[,1])-min(BCE1r[,1])),y[,1],pch=2,col=2, type = "l")
points(min(BLAr[,1])+xBLAr1*(max(BLAr[,1])-min(BLAr[,1])),y[,2],pch=3,col=3,type = "l")
points(min(BCE5r[,1])+xBCE5r1*(max(BCE5r[,1])-min(BCE5r[,1])) ,y[,3],pch=4,col=4,type = "l")
legend("topright", c("CeCO_3_F","LaCO_3_F","to be inferred"),col=c(2,3,4), lty=c(2,2,1),
       cex = 0.5)
title("(a) Experimental data")

## plot the estimated derivatives and their confidence bands.
plot(800 + x.grid * 200, y.hat[,1], type="n", xlim=c(800, 1080), ylim=c(-25, 25),
     xlab="Raman Shift", ylab="Estimated first derivatives", cex.axis=0.9)
points(800 + x.grid * 200, (2/200) * y.hat[,1], pch=2, col=2, type="l",lwd=1)
# Plot and fill the confidence band for the first series (CeCO_3_F)
#lines(800 + x.grid * 200, (2/200) * CI.BCE1$CI.lower, type="l", col=2, lty=2)
#lines(800 + x.grid * 200, (2/200) * CI.BCE1$CI.upper, type="l", col=2, lty=2)
polygon(c(800 + x.grid * 200, rev(800 + x.grid * 200)),
        c((2/200) * CI.BCE1$CI.lower, rev((2/200) * CI.BCE1$CI.upper)),
        col=rgb(1, 0, 0, alpha=0.5), border=FALSE)
points(800 + x.grid * 200, (2/200) * y.hat[,2], pch=3, col=3, type="l",lwd=1)
#lines(800 + x.grid * 200, (2/200) * CI.BLA$CI.lower, type="l", col=3, lty=2)
#lines(800 + x.grid * 200, (2/200) * CI.BLA$CI.upper, type="l", col=3, lty=2)
polygon(c(800 + x.grid * 200, rev(800 + x.grid * 200)),
        c((2/200) * CI.BLA$CI.lower, rev((2/200) * CI.BLA$CI.upper)),
        col=rgb(0, 1, 0, alpha=0.5), border=FALSE)

points(800 + x.grid * 200, (2/200) * y.hat[,3], pch=4, col=4, type="l",lwd=1)
#lines(800 + x.grid * 200, (2/200) * CI.BCE5$CI.lower, type="l", col=4, lty=2)
#lines(800 + x.grid * 200, (2/200) * CI.BCE5$CI.upper, type="l", col=4, lty=2)
polygon(c(800 + x.grid * 200, rev(800 + x.grid * 200)),
        c((2/200) * CI.BCE5$CI.lower, rev((2/200) * CI.BCE5$CI.upper)),
        col=rgb(0, 0, 1, alpha=0.5), border=FALSE)
legend("topright", c("CeCO_3_F","LaCO_3_F","to be inferred"), col=c(2,3,4),
       lty=c(2,2,1), box.lwd=1, cex=0.65)
title("(b) Est. derivatives-Iterated plugin")









par(mfrow=c(3,1))
# For CI.BCE1
plot(x.grid, CI.BCE1$Estimate, type="l")
#lines(x.grid, CI.BCE1$CI.lower, type="l", col="red", lty=2)
#lines(x.grid, CI.BCE1$CI.upper, type="l", col="red", lty=2)
polygon(c(x.grid, rev(x.grid)), c(CI.BCE1$CI.lower, rev(CI.BCE1$CI.upper)),
        col=rgb(0, 1, 0, alpha=0.4), border=F)

# For CI.BLA
plot(x.grid, CI.BLA$Estimate, type="l")
#lines(x.grid, CI.BLA$CI.lower, type="l", col="red", lty=2)
#lines(x.grid, CI.BLA$CI.upper, type="l", col="red", lty=2)
polygon(c(x.grid, rev(x.grid)), c(CI.BLA$CI.lower, rev(CI.BLA$CI.upper)),
        col=rgb(0, 1, 0, alpha=0.4), border=F)

# For CI.BCE5
plot(x.grid, CI.BCE5$Estimate, type="l")
#lines(x.grid, CI.BCE5$CI.lower, type="l", col="red", lty=2)
#lines(x.grid, CI.BCE5$CI.upper, type="l", col="red", lty=2)
polygon(c(x.grid, rev(x.grid)), c(CI.BCE5$CI.lower, rev(CI.BCE5$CI.upper)),
        col=rgb(0, 1, 0, alpha=0.4), border=F)



## Making inference
## Is it CeCO3 F?
inferredCeCO3F <- mean(abs(y.hat[,1]-y.hat[,3]))
## Is it LaCO3F?
inferredCeCO3F
inferredLaCO3F <- mean(abs(y.hat[,2]-y.hat[,3]))
inferredLaCO3F
decision<-ifelse(inferredCeCO3F< inferredLaCO3F, print("The Raman Spectrum is that of CeCO_3F"),
                 print("The Raman Spectrum is that of LaCO_3F"))

decision






## mixed model
#library(HRW); data(WarsawApts)

library(nlme)
B <- bbase(x=xBCE1r1 ,xl=min(xBCE1r1), xr=max(xBCE1r1), ndx = ndx, deg = deg)
Z=B$B
dim(Z)
y<-BCE1r[,2]
length(y)
x=xBCE1r1
length(x)
dummyID <- factor(rep(1,length(x)))

fit <- lme(y ~0,random = list(dummyID = pdIdent(~Z)))
summary(fit)
betaHat <- fit$coef$fixed ; uHat <- unlist(fit$coef$random) 
sigsqepsHat <- fit$sigma^2
sigsquHat <- as.numeric(VarCorr(fit)[1,1])
fit

dummyID <- factor(rep(1,length(B.grid$B[1])))
new_data <- data.frame(Z =B.grid$B , dummyID = dummyID)

predictions <- predict(fit, newdata =new_data  )

