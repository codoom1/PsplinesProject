# ============================================================================
# File: ComparisonAnalysis.R
# Purpose: Compare different methods for estimating derivatives of regression functions
# Author: [Christopher Odoom]
# Date: May 1, 2025
#
# Description:
#   - Defines two functions (a and b) and their derivatives.
#   - Simulates data with varying noise levels.
#   - Compares multiple derivative estimation methods (naive, plugin, resubstitution, oracle, GCPDAvg, Dai et al.).
#   - Plots mean regression functions and their derivatives.
#   - Saves simulation results for further analysis.
#
# Usage:
#   - Source this script in R or run interactively.
#   - Requires: RegSplineFunctions.R, other_methods.R, and standard R packages.
#
# Sections:
#   1. Setup and function definitions
#   2. Plotting mean regression functions and derivatives
#   3. Simulation for different noise levels and estimation methods
#   4. Saving results
# ============================================================================


# ---- Load libraries and source relevant function files ----
library(here)
library(sfsmisc)
library(locpol)
library(ggplot2)
library(dplyr)
library(tidyr)

source("R/MainFunctions.R")
source("R/competingMethods.R")

# ---- Set up parameters ----
set.seed(200)
n <- 600
x <- seq(0, 1, length.out=n)
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
ra0 <- range(fa)[2]-range(fa)[1]
ra0
x <- seq(-1, 1, length.out=n)
fb <-(sin(2*pi*x))+cos(2*pi*x)+log(4/3+x)
rb0 <- range(fb)[2]-range(fb)[1]
rb0
sf<- ra0/rb0

# ---- Plot mean regression functions and their derivatives ----
set.seed(200)
n <- 600
x <- seq(0, 1, length.out=n)
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
ra <- range(fa)[2]-range(fa)[1]
ra
fa.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
fa.pprime  <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)

x <- seq(-1, 1, length.out=n)
fb <-sf*((sin(2*pi*x))+cos(2*pi*x)+log(4/3+x) )
rb <- range(fb)[2]-range(fb)[1]
rb
fb.prime <- sf*(-2*pi*sin(2 *pi*x) + 2 *pi*cos(2*pi*x) + 1/(x + 4 / 3) )
fb.pprime <- sf*(-4*pi^2*sin(2*pi*x) - 4*pi^2*cos(2*pi*x) - 1/(x + 4 / 3)^2)

# ---- Create a white background plot ----
par(mfrow=c(2,3))
plot(x, fa, type = "l", col = "blue", xlab = "x", ylab = paste0("f.", 'a','(x)'))
plot(x, fa.prime, type="l", col = "blue", ylab = paste0("f '.", 'a','(x)'))
plot(x, fa.pprime, type="l", col = "blue", ylab = paste0("f ''.", 'a','(x)'))

plot(x, fb, col = "red", type="l", ylab = paste0("f.", 'b','(x)'))
plot(x, fb.prime, col = "red", type="l", ylab = paste0("f '.", 'b','(x)'))
plot(x, fb.pprime, col = "red", type="l", ylab = paste0("f ''.", 'b','(x)'))




# ---- Create Biometrika-style grayscale plots (new function) ----

  graphics.off()
  
  # Set up the data (same as above)
  set.seed(200)
  n <- 600
  x_a <- seq(0, 1, length.out=n)
  fa <- 32 * exp(-8 * (1 - 2 * x_a)^2) * (1 - 2 * x_a)
  fa.prime <- (4096 * x_a^2 - 4096 * x_a + 960) * exp(-8 * (1 - 2 * x_a)^2)
  fa.pprime  <- -(262144*x_a^3-393216*x_a^2+184320*x_a-26624)*exp(-8*(1-2*x_a)^2)
  
  x_b <- seq(-1, 1, length.out=n)
  fb <-sf*((sin(2*pi*x_b))+cos(2*pi*x_b)+log(4/3+x_b) )
  fb.prime <- sf*(-2*pi*sin(2 *pi*x_b) + 2 *pi*cos(2*pi*x_b) + 1/(x_b + 4 / 3) )
  fb.pprime <- sf*(-4*pi^2*sin(2*pi*x_b) - 4*pi^2*cos(2*pi*x_b) - 1/(x_b + 4 / 3)^2)
  
  
  
  # Also save as PDF
  pdf(here::here("output", "plots", "compareMethods_plots","functions_biometrika.pdf"), width = 10, height = 5)
  # Increase outer margins and space between plots
  par(mfrow=c(2,3), mar=c(5, 5, 4, 2) + 0.1, oma=c(1, 5, 0, 0), mgp=c(3, 1, 0))
  # Set up for L-shaped boxes
  par(bty = "l")
  
  # Create layout with added space for row labels
  layout(matrix(c(
    0, 1, 2, 3,
    0, 4, 5, 6
  ), nrow=2, byrow=TRUE), widths=c(0.1, 1, 1, 1))
  

  
  # Function a plots (top row) - solid lines
  plot(x_a, fa, type = "l", col = "black", lwd = 4, xlab = "x", cex.lab=1.2, 
       ylab = expression(f[a]^{(0)}*"("*x*")"), bty = "l")
  plot(x_a, fa.prime, type="l", col = "black", lwd = 4, xlab = "x", cex.lab=1.2,
       ylab = expression(f[a]^{(1)}*"("*x*")"), bty = "l")
  plot(x_a, fa.pprime, type="l", col = "black", lwd = 4, xlab = "x", cex.lab=1.2,
      ylab = expression(f[a]^{(2)}*"("*x*")"), bty = "l")
    # Add row labels before plotting
  mtext("(a)", side=2, line=.1, at=0.98, las=1, cex=2.5, outer=TRUE)

  # Function b plots (bottom row) - dashed lines
  plot(x_b, fb, type="l", col = "black", lwd = 4, lty = 2, xlab = "x", cex.lab=1.2,
       ylab = expression(f[b]^{(0)}*"("*x*")"), bty = "l")
  plot(x_b, fb.prime, type="l", col = "black", lwd = 4, lty = 2, xlab = "x", cex.lab=1.2,
       ylab = expression(f[b]^{(1)}*"("*x*")"), bty = "l")
  plot(x_b, fb.pprime, type="l", col = "black", lwd = 4, lty = 2, xlab = "x", cex.lab=1.2,
       ylab = expression(f[b]^{(2)}*"("*x*")"), bty = "l")
  mtext("(b)", side=2, line=.1, at=0.45, las=1, cex=2.5, outer=TRUE)
  dev.off()
  
  cat("Biometrika-style function plots saved as:\n")
  cat("- functions_biometrika.pdf\n")



# ---- Simulation for different noise levels and function a ----
noise.sim <- function(sig= 0.01, n=300, MC.N=10, r=1,bdeg=3,pord=2, 
                        funtype=c("a","b"),hs=seq(0.01,0.99,0.05)){
  ## define the function 
  if(funtype=="a"){
    x <- seq(from=0,to=1,length=n)
    x.grid <- x
    f <-32*exp(-8*(1-2*x)^2)*(1-2*x)
    if(r==0){
      fr.grid <- 32*exp(-8*(1-2*x.grid)^2)*(1-2*x.grid)
    } else if(r==1){
      fr.grid<-  (4096*x.grid^2-4096*x.grid+960)*exp(-8*(1-2*x.grid)^2)
    }else{
      fr.grid<- -(262144*x.grid^3-393216*x.grid^2+184320*x.grid-26624)*exp(-8*(1-2*x.grid)^2)
    }
    plot(x, f, type = "l", col="red")
    
  }else if(funtype=="b"){
    x <- seq(from=-1,to=1,length=n)
    x.grid <- x
    f<- sf*(sin(2*pi*x)+cos(2*pi*x)+log(4/3+x))
    if(r==0){
      fr.grid <- sf*(sin(2*pi*x.grid)+cos(2*pi*x.grid)+log(4/3+x.grid))
    } else if(r==1){
      fr.grid<-  sf*(-2*pi*sin(2 *pi*x.grid) + 2 *pi*cos(2*pi*x.grid) + 1/(x.grid + 4 / 3))
    }else{
      fr.grid<- sf*(-4*(pi^2)*sin(2*pi*x.grid) - 4*(pi^2)*cos(2*pi*x.grid) - 1/(x.grid + 4 / 3)^2)
    }
    plot(x, f, type = "l", col="red")
  }
  
  sig<- sig ## amount of noise
  r <- r
  nseg <- 35
  bdeg <- bdeg
  
  keep.emise <- data.frame(
                        naive=rep(NA,MC.N),
                        simple =rep(NA,MC.N),
                        plug.in=rep(NA,MC.N),
                        resub=rep(NA,MC.N),
                        GCPDAvg=rep(NA,MC.N),
                        Dailc=rep(NA, MC.N),
                        oracle=rep(NA,MC.N))
  
  
  for (i in 1:MC.N)
  {	
    
    y <- f + sig*rnorm(n)
    naive.fit <- naive.est.opt(x = x, y = y, r = r, x.grid = x.grid, nseg = nseg, pord = pord, bdeg = bdeg) # nolint

    simple.fit <- simple.est(lambda=naive.fit$lambda , r=r, x=x,x.grid, f= naive.fit$f.hat,bdeg=bdeg, nseg=nseg,pord=pord) # nolint

    plugin.fit <- plugin.est(x=x, y=y, r=r,nseg=nseg , pord = pord ,
                        bdeg=bdeg, x.grid=x.grid,fr="Y") # nolint

    resub.fit <- resub.est(x = x, y = y, r = r, x.grid = x.grid, nseg= nseg, pord = pord,bdeg = bdeg, tol = 1e-5, ITs = 100) # nolint

    oracle.fit <- oracle.est(initial.lambda = resub.fit$lambda, x = x, y = y, r = r, fr.grid = fr.grid, nseg= nseg, pord = pord, bdeg = bdeg, x.grid = x.grid) # nolint

    GCPDAvg <- charnigo.est(x=x, y=y, r=r, x.grid=x.grid,k=c(5,10,15,20,25),k1=c(5,15,25),
                                        k2=c(5,15,25),nknots=nseg, hs=hs)
    if(r==1){
    Daifit<-  dai.fit(x, y,p=r,q=7,x.grid, fr.grid=fr.grid )
    }else{
      Daifit<-  dai.fit(x, y,p=r,q=8,x.grid, fr.grid=fr.grid )
     
    }
    
    
    keep.emise$naive[i]<-  mean((fr.grid-naive.fit$fr.hat )^2)
    keep.emise$simple[i]<-  mean((fr.grid-simple.fit$fr.hat )^2)
    keep.emise$plug.in[i]<-mean((fr.grid-plugin.fit$fr.hat )^2)
    keep.emise$resub[i]<-  mean((fr.grid-resub.fit$fr.hat )^2)
    keep.emise$oracle[i]<-  mean((fr.grid-oracle.fit$fr.hat)^2)
    
    keep.emise$GCPDAvg[i] <- mean((GCPDAvg$fr.hat  -fr.grid)^2)
    keep.emise$Dailc[i] <- mean((Daifit$fr.hat -fr.grid)^2)
    print(GCPDAvg$h)
    if (i %% 10 == 0) {

    png(filename = here::here("output", "plots", "compareMethods_plots",paste0("compare",i,".png")), width = 800, height = 600)
    plot(x.grid, fr.grid, type = "l", col="black", main = "Deriv.Est")
    lines(x.grid, naive.fit$fr.hat, col="blue", main="naive")
    lines(x.grid, simple.fit$fr.hat, col="#a2a3a3", main="simple")
    lines(x.grid, plugin.fit$fr.hat, col="red", main="plugin")
    lines(x.grid, resub.fit$fr.hat, col="green", main="resub")
    lines(x.grid, oracle.fit$fr.hat, col="wheat", main="oracle")
    lines(x.grid, GCPDAvg$fr.hat, type = "l", col = "gold")
    lines(x.grid, Daifit$fr.hat, type = "l", col = "darkblue")
    legend("bottomright", legend = c("true", "naive","simple", "plugin", "resub","oracle", "GCPDAvg",
                                     "Dai_etal"),
           col = c("black", "blue","#a2a3a3", "red", "green", "wheat","gold","darkblue"), lty = 1,box.lwd=0.2,cex=0.3)
    dev.off()
    }
    print("*******")
    print(keep.emise[i,])
    print("******")
    
  }
  print("*******")
  print(sig)
  print("*******")
  print(keep.emise)
  return(keep.emise)
} 

# ---- Run all simulations and save results ----
run_all_simulations <- function(sigs = c(0.1, 0.5, 2), MC.N = 10, n = 100) {
  # First derivative
  r <- 1
  bdeg <- 4
  pord <- 2
  simb1 <- lapply(sigs, function(sig) {
    noise.sim(sig = sig, n = n, MC.N = MC.N, r = r, bdeg = bdeg, pord = pord,
              funtype = "b", hs = seq(.01, 0.5, length.out = 21))
  })
  save(simb1, file = here::here("output", "data", "compareMethods_data", "emseb1.RData"))
  print("Finished: First derivative, function b (emseb1.RData)")

  sima1 <- lapply(sigs, function(sig) {
    noise.sim(sig = sig, n = n, MC.N = MC.N, r = r, bdeg = bdeg, pord = pord,
              funtype = "a", hs = seq(.01, 0.5, length.out = 21))
  })
  save(sima1, file = here::here("output", "data", "compareMethods_data", "emsea1.RData"))
  print("Finished: First derivative, function a (emsea1.RData)")

  # Second derivative
  r <- 2
  bdeg <- 4
  pord <- 3
  simb2 <- lapply(sigs, function(sig) {
    noise.sim(sig = sig, n = n, MC.N = MC.N, r = r, bdeg = bdeg, pord = pord,
              funtype = "b", hs = seq(.01, 0.5, length.out = 21))
  })
  save(simb2, file = here::here("output", "data", "compareMethods_data", "emseb2.RData"))
  print("Finished: Second derivative, function b (emseb2.RData)")

  sima2 <- lapply(sigs, function(sig) {
    noise.sim(sig = sig, n = n, MC.N = MC.N, r = r, bdeg = bdeg, pord = pord,
              funtype = "a", hs = seq(.01, 0.5, length.out = 21))
  })
  save(sima2, file = here::here("output", "data", "compareMethods_data", "emsea2.RData"))
  print("Finished: Second derivative, function a (emsea2.RData)")
}

# ---- Run all simulations with specified parameters ----
run_all_simulations(MC.N = 600, n = 600)


## End of simulation script------##



