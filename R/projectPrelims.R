# This script performs penalized spline analysis and visualization of regression functions and their derivatives.
# It includes the following functionalities:
# - Loading necessary libraries and external functions.
# - Generating regression functions and their derivatives.
# - Visualizing regression functions and their derivatives.
# - Implementing and testing Generalized Cross-Validation (GCV) for selecting smoothing parameters.
# - Comparing different estimation approaches for Mean Integrated Squared Error (MISE).
# - Investigating the behavior of smoothing matrices and their properties.
# - Conducting Monte Carlo simulations to evaluate estimation methods.
# - Exploring the effect of sample size and smoothing parameters on estimation accuracy.
# - Utilizing penalized splines as kernel smoothers and analyzing their behavior.

## using the functions
#rm(list = ls())


## Load the required libraries and functions
library(here)
library(dplyr)  # For data manipulation
library(ggplot2)  # For data visualization
library(MASS)  # For statistical functions
library(lattice)  # For advanced plotting
library(gridExtra)  # For arranging multiple plots
library(mgcv)  # For generalized additive models
library(tidyr)  # For data tidying

# Source external functions from "RegSplineFunctions.R"
source(here::here("R", "mainFunctions.R"))

# Set seed for reproducibility
set.seed(200)

# Define the number of data points
n <- 600

# Generate the first regression function (fa) and calculate its range
x <- seq(0, 1, length.out=n)  # Generate a sequence of x values
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)  # Define the function fa
ra0 <- range(fa)[2] - range(fa)[1]  # Calculate the range of fa
ra0  # Print the range of fa

# Generate the second regression function (fb) and calculate its range
x <- seq(-1, 1, length.out=n)  # Generate a sequence of x values
fb <- (sin(2 * pi * x)) + cos(2 * pi * x) + log(4 / 3 + x)  # Define the function fb
rb0 <- range(fb)[2] - range(fb)[1]  # Calculate the range of fb
rb0  # Print the range of fb

# Calculate the scaling factor to normalize the functions
sf <- ra0 / rb0

#### Making plots of the mean regression functions and their derivatives
# Generate the first regression function and its derivatives
x <- seq(0, 1, length.out=n)  # Generate a sequence of x values
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)  # Define the function fa
ra <- range(fa)[2] - range(fa)[1]  # Calculate the range of fa
ra  # Print the range of fa
fa.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)  # First derivative of fa
fa.pprime <- -(262144 * x^3 - 393216 * x^2 + 184320 * x - 26624) * exp(-8 * (1 - 2 * x)^2)  # Second derivative of fa

# Generate the second regression function and its derivatives
x <- seq(-1, 1, length.out=n)  # Generate a sequence of x values
fb <- sf * ((sin(2 * pi * x)) + cos(2 * pi * x) + log(4 / 3 + x))  # Define the function fb
rb <- range(fb)[2] - range(fb)[1]  # Calculate the range of fb
rb  # Print the range of fb
fb.prime <- sf * (-2 * pi * sin(2 * pi * x) + 2 * pi * cos(2 * pi * x) + 1 / (x + 4 / 3))  # First derivative of fb
fb.pprime <- sf * (-4 * pi^2 * sin(2 * pi * x) - 4 * pi^2 * cos(2 * pi * x) - 1 / (x + 4 / 3)^2)  # Second derivative of fb

# Create a white background plot
# Save the plot as a PNG file
png(here::here("output", "plots","prelims_plots", "functions.png"), width = 800, height = 600)
par(mfrow=c(2,3))
x <- seq(0, 1, length.out=n)
plot(x, fa, type = "l", col = "blue", xlab = "x", ylab = paste0("f.", 'a','(x)'))
plot(x, fa.prime, type="l", col = "blue", ylab = paste0("f '.", 'a','(x)'))
plot(x, fa.pprime, type="l", col = "blue", ylab = paste0("f ''.", 'a','(x)'))

x <- seq(-1, 1, length.out=n)
plot(x, fb, col = "red", type="l", ylab = paste0("f.", 'b','(x)'))
plot(x, fb.prime, col = "red", type="l", ylab = paste0("f '.", 'b','(x)'))
plot(x, fb.pprime, col = "red", type="l", ylab = paste0("f ''.", 'b','(x)'))

dev.off()
par(mfrow=c(1,1)) 

### Testing core functions in "RegSplineFunctions.R"

########### Testing the GCV criterion function ############
### Simulate data for testing
#rm(list = ls())
set.seed(200)
n <- 600
x <- seq(0, 1, length.out=n)
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
fa.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
fa.pprime  <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)

sig<- .5 ## noise
f<-fa
y <- f + sig*rnorm(n) ## model  fx +epsilon
r=1
nseg = 35
pord = 2
bdeg = 4

##########
## Using GCV to select a \lambda to est f.

l.est<- gcvlambda(lambda = 1, x=x,y=y,
                   nseg=35, pord = 2, bdeg = 4)

log.lam <- seq(-2, 2, length.out=500)

lm1<- 10^log.lam
print(lm1)

searh<- sapply(lm1, function(l){
  gcvlambda(lambda = l, x=x,y=y,
            nseg=nseg, pord = pord, bdeg = bdeg)
  
})
idx<- which.min(searh) ## find the index of the minimum

# Save the plot as a PNG file
png(here::here("output", "plots","prelims_plots", "gcv.png"), width = 800, height = 600)
plot(lm1, (searh), type="l", col="blue", lwd=2,
     xlab = "log10(lambda)", ylab = "log10GCV")

dev.off()
################### End of GCV #######################

### This function compares the different type of objective
## functions
lams=lm1
trial <- sig.effect(sig=.09, lams=lm1, x=x,r=1, f=fa, fr=fa.prime)
mise.true<- trial$mise.true
mise.est_ideal <- trial$mise.est_ideal
mise.est_prop<- trial$mise.est_prop

# Save the plot as a PNG file
png(here::here("output", "plots","prelims_plots", "plot_VarSqbias.png"), width = 800, height = 600)
par(mfrow = c(1, 3))
# First plot
plot((lams), (unlist(mise.true["var",])), col = "black", type = "l",
     xlab = "lambda", ylab = "var", ylim = c(0, 2.5))
lines((lams), (unlist(mise.est_ideal["var",])), col = "red", type = "l")
lines((lams), (unlist(mise.est_prop["var",])), col = "blue", type = "l")
legend("topright", legend = c("True Var", "original Est", "Propose Est"),
       col = c("black", "red", "blue"), lty = 1, cex = 1)

# Second plot
plot((lams), (unlist(mise.true["sq.bias",])), col = "black", type = "l",
     xlab = "lambda", ylab = "sq.bias", ylim = c(0, 17))
lines((lams), (unlist(mise.est_ideal["sq.bias",])), col = "red", type = "l")
lines((lams), (unlist(mise.est_prop["sq.bias",])), col = "blue", type = "l")
legend("topleft", legend = c("True Sq.Bias", "original Est", "Propose Est"),
       col = c("black", "red", "blue"), lty = 1, cex = 1)

# Third plot
plot((lams), (unlist(mise.true["mise",])), col = "black", type = "l",
     xlab = "lambda", ylab = "mise", ylim = c(0, 17))
lines((lams), (unlist(mise.est_ideal["mise",])), col = "red", type = "l")
lines((lams), (unlist(mise.est_prop["mise",])), col = "blue", type = "l")
legend("topleft", legend = c("True mise", "original Est", "Propose Est"),
       col = c("black", "red", "blue"), lty = 1, cex = 1)

dev.off()  # Close the device

# Reset layout to default
par(mfrow = c(1, 1))

ind11<- which.min((unlist(mise.est_prop["mise",])))
print(paste("Optimal lambda for proposed estimator:", lams[ind11]))
ind12<- which.min((unlist(mise.est_ideal["mise",])))
print(paste("Optimal lambda for original estimator:", lams[ind12]))

######### End of Comparison of approaches for the MISE ########

############# Is lambda resub same as minimizer of mise.lambda #############
### plugin estimate of fr

pin.esty <-  plugin.estfr(x=x, y=y, r=r,nseg=nseg, pord = pord,
                        bdeg=bdeg, x.grid=x,fr="Y")
pin.esty$lambda
pin.estn <-  plugin.estfr(x=x, y=y, r=r,nseg=nseg, pord = pord,
                          bdeg=bdeg, x.grid=x,fr=NULL)

pin.estn$lambda

### Using the naive estimator to get a derivative est and f.hat
n.est <- naive.est.opt(x=x, y=y,r=r
                   , nseg = nseg, bdeg = bdeg,pord = pord,
                   x.grid = x)
n.est$fr.est$lambda ### Resubstitution estimate of fr.
resub.fr<- resub.est(x=x, y=y, r=r, x.grid = x,nseg=nseg, pord = pord,
                     bdeg=bdeg,tol=1e-10, ITs=100)
resub.fr$lambda ## The output here shows that the resub converges

## estimate of the mean function for a given lambda
png(here::here("output", "plots","prelims_plots", "plot_estimates.png"), width = 800, height = 600)
par(mfrow=c(1,2))
plot(x, f, type = "l", col="red")
lines(x, n.est$f.hat, type = "l")
lines(x, pin.esty$f.hat, type = "l", col="blue")
lines(x, pin.estn$f.hat, type = "l", col="yellow")
legend("topright", legend = c("True f", "Naive est", "Plugin estY", "Plugin estN"),
       col = c("red", "black", "blue", "yellow"), lty = 1, cex = 1)

plot(x, fa.prime,type = "l", col="red")
lines(x, n.est$fr.hat, type = "l")
lines(x, pin.esty$fr.hat, type = "l", col="blue")
lines(x, pin.estn$fr.hat, type = "l", col="yellow")
lines(x, resub.fr$fr.hat, type = "l", col="green")
legend("topright", legend = c("True fr", "Naive est", "Plugin estY","Plugin estN", "Resub est"),
       col = c("red", "black", "blue", "yellow", "green"), lty = 1, cex = 1)

dev.off()
par(mfrow=c(1,1)) 

##### hatmatrix behaves like a kernel #############
H<- n.est$fr.est$A
dim(H)
K <- n.est$fr.est$K
dim(K)

### Thinking of the penalized splines as kernels smoothing ####
png(here::here("output", "plots","prelims_plots", "plot_kernels.png"), width = 800, height = 600)
par(mfrow=c(1,2))
plot(1:600, H[300,], type = "l")
plot(1:600, K[300,], type = "l", col="red")
legend("topright", legend = c("H", "K"),
       col = c("black", "red"), lty = 1, cex = 1)
dev.off()

############# Compare the methods based on lambda chosen ############

compare<- sig.lambda.fun(sig=0.08,f=fa,fr=fa.prime,x=x, x.grid=x, r=r, nseg=nseg, pord=pord, 
                  bdeg=bdeg ,tol=1e-10, ITs=100)

compare

############ End of Comparison ###################

######### Checking  the  effect of lambda on Kernels #########
emise.1 <- mise.lambda(lambda=.1, x=x,y=y, r=1,sig=0.1,nseg=35, pord=2, 
                        bdeg=4, f=f, fr=NULL)
emise100 <- mise.lambda(lambda=100, x=x,y=y, r=1,sig=0.1,nseg=35, pord=2, 
                     bdeg=4, f=f, fr=NULL)
emise10000000 <- mise.lambda(lambda=10000000, x=x,y=y, r=1,sig=0.1,nseg=35, pord=2, 
                     bdeg=4, f=f, fr=NULL)

print(paste("lambda = 0.1, sum of abs diag H:", sum(abs(diag(emise.1$H)))))
png(here::here("output", "plots","prelims_plots", "plot_kernels2.png"), width = 800, height = 600)
par(mfrow=c(1,3))
plot(1:600, emise.1$H[300,], type = "l", col="red")
plot(1:600, emise100$H[300,], type = "l", col="red")
plot(1:600, emise10000000$H[300,], type = "l", col="red")
dev.off()
par(mfrow=c(1,1))

############ End of Checking Kernel ###############

######### Checking the effect of lambda on the MISE #########
ls<- 1:100
store <- numeric()
mise.null<- numeric()
mise<- numeric()
for(i in 1:length(ls) ){
  
emise.1 <- mise.lambda(lambda=ls[i], x=x,y=y, r=1,sig=0.5,
                       nseg=35, pord=2, 
                       bdeg=4, f=f, fr=NULL)
emise.2 <- mise.lambda(lambda=ls[i], x=x,y=y, r=1,sig=0.5,
                       nseg=35, pord=2, 
                       bdeg=4, f=f, fr=fa.prime)
store[i]<- sum(abs(diag(emise.1$H)))
mise.null[i] <- emise.1$mise
mise[i] <- emise.2$mise
}

png(here::here("output", "plots","prelims_plots", "plot_mise_compare.png"), width = 800, height = 600)
par(mfrow=c(1,2))
plot(ls, store, type = "l", main = "abs.sum H")
plot(ls, mise.null, type="l",main = " comparison of two ests")
lines(ls, mise, col="red")
dev.off()
par(mfrow=c(1,1))
store[99:100]
mise.null[99:100]
mise[99:100]
emise.1$mise
emise.1$var
emise.1$sq.bias

######### simple estimate of the derivative ########
### Hf
##test function ###
set.seed(299)
n <- 900
x <- seq(0, 1, length.out=n)
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
fa.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
fa.pprime  <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
x.grid <- seq(0, 1, length.out=300)
sig<- .01 ## noise
f<-fa
fr <- fa.prime
y <- f + sig*rnorm(n) ## model  fx +epsilon
r=1
nseg = 37
pord = 1
bdeg = 3

n.est <- naive.est.opt(x=x, y=y,r=r
                   , nseg = nseg, bdeg = bdeg,pord = pord,
                   x.grid = x)
n.est$fr.est$lambda

est.fr <- simple.est(lambda=n.est$fr.est$lambda , r=r, x=x,x.grid=x.grid, 
             f= n.est$fr.est$f.hat,bdeg=bdeg, nseg=nseg)

png(here::here("output", "plots", "prelims_plots", "plot_estimates2.png"), width = 800, height = 600)
par(mfrow=c(1,2))
plot(x, fr, type = "l", col="blue")
lines(x,est.fr$simp.estfr, col="red" )
lines(x, n.est$fr.est$fr.hat)

plot(x, n.est$fr.est$fr.hat, type="l")
lines(x,est.fr$simp.estfr, col="red" )
dev.off()
par(mfrow=c(1,1))
mean((fr-est.fr$simp.estfr)^2)
mean((fr-n.est$fr.est$fr.hat)^2)

#### MC samples to compare simple estimate of the derivative and the naive estimate #######
set.seed(20)
n <- 600
x <- seq(0, 1, length.out=n)
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
fa.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
fa.pprime  <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
frs<- matrix(c(fa.prime, fa.pprime), ncol = 2)
x.grid <- seq(0, 1, length.out=300)
sig<- .1 ## noise level
f<-fa ## true function
nsim<- 10

keep.mise <- data.frame(
  r = rep(NA, nsim*2),
  mise.simp.opt = rep(NA, nsim*2),
  mise.simp.grid = rep(NA, nsim*2),
  mise.naive.opt = rep(NA, nsim*2),
  mise.naive.grid = rep(NA, nsim*2)
)
counter<- 1
for(r in 1:2){
  fr <- frs[,r]
  for(i in 1:nsim){
    y <- f + sig*rnorm(n)
    n.est <- naive.est.opt(x=x, y=y,r=r
                           , nseg = nseg, bdeg = bdeg,pord = pord,
                           x.grid = x)
    n.est1 <- naive.est(x=x, y=y,r=r,log.lam= seq(-2.5, 2, length.out=100)
                        , nseg = nseg, bdeg = bdeg,pord = pord,
                        x.grid = x)
    est.fr <- simple.est(lambda=n.est$lambda , r=r, x=x,x.grid=x.grid, 
                 f= n.est$f.hat,bdeg=bdeg, nseg=nseg)
    
    est.fr1 <- simple.est(lambda=n.est1$fr.est$lambda , r=r, x=x,x.grid=x.grid, 
                  f= n.est1$fr.est$f.hat,bdeg=bdeg, nseg=nseg)
    
    keep.mise$mise.simp.opt[counter] <-  mean((fr-est.fr$simp.estfr)^2)
    keep.mise$mise.simp.grid[counter] <-  mean((fr-est.fr1$simp.estfr)^2)
    keep.mise$mise.naive.opt[counter] <- mean((fr-n.est$fr.hat)^2)
    keep.mise$mise.naive.grid[counter] <- mean((fr-n.est1$fr.est$fr.hat)^2)
    keep.mise$r[counter] <- r
    counter<- counter+1
    print(i)
  }
  print("*****")
  print(r)
}

keep.misehead <- head(keep.mise)
print(keep.misehead)
# Save the plot as a PNG file
keep.mise.long <- keep.mise %>%
  pivot_longer(cols = c(mise.simp.opt,mise.simp.grid, mise.naive.opt,
                        mise.naive.grid), 
               names_to = "emise_type", 
               values_to = "emise_value")

# Create a boxplot for r = 1
png(here::here("output", "plots","prelims_plots", "plot_miseC.png"), width = 800, height = 600)
p1 <- ggplot(subset(keep.mise.long, r == 1), aes(x = emise_type, y = emise_value, fill = emise_type)) +
  geom_boxplot() +
  labs(
    x = "MISE Type", 
    y = "EMISE Value", 
    title = "EMISE Values for r = 1"
  ) +
  theme_bw() +
  theme(legend.position = "none")  # Remove legend for individual plots

# Create a boxplot for r = 2
p2 <- ggplot(subset(keep.mise.long, r == 2), aes(x = emise_type, y = emise_value, fill = emise_type)) +
  geom_boxplot() +
  labs(
    x = "MISE Type", 
    y = "EMISE Value", 
    title = "EMISE Values for r = 2"
  ) +
  theme_bw() +
  theme(legend.position = "none")  # Remove legend for individual plots

# Arrange the two plots side by side
grid.arrange(p1, p2, ncol = 2)
dev.off()

##### plot of estimate of the expected difference squared of 
# the true and est mise  for a fixed lambda#####
log.n<- seq(log10(1000), log10(10000),log10(1.11) )
l8<- exp.mise(lambda=.2,r=1, bdeg=4, pord=2,sig=.1,nsim=1, log.n=log.n)
emise.diff
l9<- emise.diff(lambda=2,r=1, bdeg=4, pord=2,sig=.1,nsim=1, log.n=log.n)

save(l5, file = here::here("output", "data","prelims_data", "oldrate3.RData"))
l5<- get(load(here::here("output", "data","prelims_data", "oldrate3.RData")))
floor(ns/20)
keep.difmise=l8$keep.difmise1

lmodel <- lm(log10(exp.biasdiff)~log10(n), data = keep.difmise)
summary(lmodel)
confint(lmodel)

## var difference
lmodel1 <- lm(log10(exp.vardiff)~log10(n), data = keep.difmise)
summary(lmodel1)
confint(lmodel1)
## Sigma^2 rate
lmodel2 <- lm(log10(sig.hats.diff)~log10(n), data = keep.difmise)
summary(lmodel2)
confint(lmodel2)

lmodel3 <- lm(log10(hh.sqp)~log10(n), data = keep.difmise)
summary(lmodel3)
confint(lmodel3)

#save(keep.difmise, file = "diff.RData")
#keep.difmise<- get(load("diff.RData"))

png(here::here("output", "plots","prelims_plots", "plot_mise_diff.png"), width = 800, height = 600)
p1<- keep.difmise %>% 
  ggplot(aes(x = log10(n), y = log10(exp.biasdiff))) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add regression line
  theme_bw()
p1
p2<- keep.difmise %>% 
  ggplot(aes(x = log10(n), y = log10(exp.vardiff))) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add regression line
  theme_bw()
p2
p3<- keep.difmise %>% 
  ggplot(aes(x = log10(n), y = log10(sig.hats.diff))) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add regression line
  theme_bw()
p3

p4<- keep.difmise %>% 
  ggplot(aes(x = log(n), y = log(hh.sq))) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, color = "blue") +  # Add regression line
  theme_bw()
p4

# p2<- keep.difmise %>% ggplot(aes(x=log10(n), y=log10(exp.sqdiff)))+
#   geom_line()+
#   theme_bw()
p1 <- p1 + xlab("log10(n)") + ylab("log10(var.sqdiff)")
p2 <- p2 + xlab("log10(n)") + ylab("log10(abs.biasdiff)")
p1 <- p1 + ggtitle("exp.Sq diff(Variance)")
p2 <- p2 + ggtitle("exp.abs diff(bias)")
g <- grid.arrange(p1, p2, ncol = 2, nrow = 1)
ggsave(here::here("output", "plots","prelims_plots", "combined_plot.png"), g, width = 10, height = 5)
dev.off()

# Keep K by taking the first value per group

keep.difmise %>%group_by(n) %>%  summarise( edif = mean(exp.sqdiff))

################## Using the GCV package #########

set.seed(20)
n <- 4000
x <- seq(0, 1, length.out=n)
fa <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
fa.prime <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
fa.pprime  <- -(262144*x^3-393216*x^2+184320*x-26624)*exp(-8*(1-2*x)^2)
frs<- matrix(c(fa.prime, fa.pprime), ncol = 2)
x.grid <- seq(0, 1, length.out=300)
sig<- .1 ## noise level
f<-fa ## true function
y<- f+rnorm(n=n, sd=sig)
r<-1
nseg<- 300
bdeg <- 3
pord <- 2

f.hat <- gam(y~s(x,k=(nseg+bdeg), bs="ps", m=c(bdeg-1, pord)),
             method="GCV.Cp")
sig_hat <- sigma(f.hat)
sig_hat
tr <- sigma.est(x=x, y=y,r=r
                , nseg = nseg, bdeg = bdeg,pord = pord,
                x.grid = x)

tr

plot(x, f.hat$fitted.values, main = "estimate of the mean function",
     ylab = "f.hat",type="l")

f.hat.naive <- naive.est.opt(x=x, y=y,r=r
                             , nseg = nseg, bdeg = bdeg,pord = pord,
                             x.grid = x)

lines(x, f.hat.naive$f.hat, col="red")
lines(x, f, col="blue")
f.hat.naive$f.hat[500]
f.hat$fitted.values[500]

abs(f.hat.naive$f.hat[200]-f.hat$fitted.values[200])
mean((f.hat.naive$f.hat-f.hat$fitted.values))

f.hat$optimizer
a<- bbase(x, nseg = nseg, bdeg = bdeg)
dim(a)

ns <- samp(range_start=500,range_end=3000,points_per_range=2)
examp <- sig.sim(r=1, bdeg=3, pord=1,sig=.1,nsim=5,ns)
examp$keep.sigma

examp$keep.sigma %>% ggplot(aes(x=K, y=trace))+
  geom_point(col="red")

examp$keep.sigma %>% ggplot(aes(x=K, y=edf))+
  geom_point(col="red")

examp$keep.sigma %>% ggplot(aes(x=log10(n), y=log(sigma.gcvdiff) ))+
  geom_point(col="red")

mdl <- lm(log(sigma.gcvdiff)~log10(n), data =examp$keep.sigma  )
summary(mdl)

################ Using loops ###########

r=2
# nseg = 37
pord = 3
bdeg = 5
sig <- 0.1
nsim <- 100
lambda <-.002

# Define the range and number of samples per thousand
start <- 1000
end <- 10000
samples_per_thousand <- 2  # Adjust this as needed

# Generate the logarithmic grid for each thousand
ns <- unlist(lapply(seq(start, end - 1000, by = 1000), function(x) {
  10^seq(log10(x), log10(x + 1000), length.out = samples_per_thousand + 1)[-(samples_per_thousand + 1)]
}))

# Round to integers (if needed)
ns <- ceiling(ns)

log.n <- log10(ns)
# # Define the sample sizes and smoothing parameters
# log.n <- seq(log10(1000), log10(10000), log10(1.11))  # Logarithmic grid for n
# ns <- ceiling(10^log.n)  # Convert log scale to actual sample sizes
v <- (2 * (2 * pord)) + 1  # Smoothness parameter
rate <- 1 / v
Ks <- ceiling(10^(1.9 + rate * log.n))  # Number of segments
Ks
keep.difmise <- data.frame(
  n = rep(NA, nsim*length(ns)),
  K = rep(NA, nsim*length(ns)),
  var.diff= rep(NA, nsim*length(ns)),
  varp.diff= rep(NA, nsim*length(ns)),
  sq.bias.diff= rep(NA, nsim*length(ns)),
  mise.diff = rep(NA, nsim*length(ns)),
  sig.sqdiff = rep(NA, nsim*length(ns)),
  hh.sq = rep(NA, nsim*length(ns)),
  hh.sqp = rep(NA, nsim*length(ns))
)
target_x <- 0.5
counter <- 1
for (j in 1:length(ns)) {
  
  x <- sort(unique(c(seq(0, 1, length.out = ns[j]), target_x)))
  f <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
  
  if (r == 1) {
    fr <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
  } else if (r == 2) {
    fr <- -(262144 * x^3 - 393216 * x^2 + 184320 * x - 26624) * exp(-8 * (1 - 2 * x)^2)
  }
  
  for (i in 1:nsim) {
  
    y <- f + rnorm(n = ns[j], mean = 0, sd = sig)
    t.emise <- mise.diff(lambda = lambda, x, y, r = r, sig = sig, nseg = Ks[j], 
                         pord = pord, bdeg = bdeg, f = f, fr = fr, target_x=target_x)
    
    keep.difmise$n[counter] <- ns[j]
    keep.difmise$K[counter] <- Ks[j]
    keep.difmise$var.diff[counter] <- t.emise$var.diff
    keep.difmise$varp.diff[counter] <- t.emise$var.diffp
    keep.difmise$sq.bias.diff[counter] <- t.emise$sq.bias.diff
    keep.difmise$mise.diff[counter] <- t.emise$mise.diff
    keep.difmise$sig.sqdiff[counter] <- t.emise$sig.sqdiff
    keep.difmise$hh.sq[counter] <- t.emise$HH
    keep.difmise$hh.sqp[counter] <- t.emise$hh
    
    counter <- counter + 1
    print(i)
  }
  print(n[j])
  # Summarizing data after each `n`
  keep.difmise1 <- keep.difmise %>% 
    group_by(n) %>% 
    summarise(
      K = first(K),
      var.diff = mean(var.diff, na.rm = TRUE),
      varp.diff = mean(varp.diff, na.rm = TRUE),
      sq.bias.diff = mean(sq.bias.diff, na.rm = TRUE),
      mise.diff = mean(mise.diff, na.rm = TRUE),
      sig.sqdiff = mean(sig.sqdiff, na.rm = TRUE),
      hh.sq = mean(hh.sq, na.rm = TRUE),
      hh.sqp = mean(hh.sqp, na.rm = TRUE)
    )
  keep.difmise1 <- keep.difmise1 %>%
    filter(!is.na(mise.diff), !is.na(var.diff),
           !is.na(varp.diff),!is.na(sq.bias.diff), !is.na(sig.sqdiff))
  # Compute slopes
  slope_mise <- coef(lm(log10(mise.diff) ~ log10(n), data = keep.difmise1))[2]
  slope_var <- coef(lm(log10(var.diff) ~ log10(n), data = keep.difmise1))[2]
  slope_varp <- coef(lm(log10(varp.diff) ~ log10(n), data = keep.difmise1))[2]
  slope_bias <- coef(lm(log10(sq.bias.diff) ~ log10(n), data = keep.difmise1))[2]
  slope_sig <- coef(lm(log10(sig.sqdiff) ~ log10(n), data = keep.difmise1))[2]
  slope_HH <- coef(lm(log10(hh.sq) ~ log10(n), data = keep.difmise1))[2]
  slope_hh <- coef(lm(log10(hh.sqp) ~ log10(n), data = keep.difmise1))[2]
  # Print slopes in console for debugging
  print(paste("Slope of mise.diff:", round(slope_mise, 3)))
  print(paste("Slope of var.diff:", round(slope_var, 3)))
  print(paste("Slope of var.diffp:", round(slope_varp, 3)))
  print(paste("Slope of sq.bias.diff:", round(slope_bias, 3)))
  print(paste("Slope of sig.sqdiff:", round(slope_sig, 3)))
  print(paste("Slope of HH:", round(slope_HH, 3)))
  print(paste("Slope of hh:", round(slope_hh, 3)))
  # Plotting with slopes
  p1 <- ggplot(keep.difmise1, aes(x = log10(n), y = log10(mise.diff))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue")+ 
    # + 
    #   annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$mise.diff)), 
    #            label = paste0("Slope: ", round(slope_mise, 2)), hjust = 0, color = "red",vjust = 1) +
    theme_bw() + xlab("log10(n)") + ylab("log10(mise.diff)")
  
  p2 <- ggplot(keep.difmise1, aes(x = log10(n), y = log10(var.diff))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
    # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$var.diff)), 
    #          label = paste0("Slope: ", round(slope_var, 2)), hjust = 0, color = "red",vjust = 1) +
    theme_bw() + xlab("log10(n)") + ylab("log10(var.diff)")
  
  p3 <- ggplot(keep.difmise1, aes(x = log10(n), y = log10(varp.diff))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
    # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$var.diff)), 
    #          label = paste0("Slope: ", round(slope_var, 2)), hjust = 0, color = "red",vjust = 1) +
    theme_bw() + xlab("log10(n)") + ylab("log10(varp.diff)")
  
  p4 <- ggplot(keep.difmise1, aes(x = log10(n), y = log10(sq.bias.diff))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
    # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$sq.bias.diff)), 
    #          label = paste0("Slope: ", round(slope_bias, 2)), hjust = 0,vjust = 1, color = "red", size = 5) +
    theme_bw() + xlab("log10(n)") + ylab("log10(sq.bias.diff)")
  
  p5 <- ggplot(keep.difmise1, aes(x = log10(n), y = log10(sig.sqdiff))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
    # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$sig.sqdiff)), 
    #          label = paste0("Slope: ", round(slope_sig, 2)), hjust = 0, color = "red",vjust = 1) +
    theme_bw() + xlab("log10(n)") + ylab("log10(sig.sqdiff)")
  
  p6 <- ggplot(keep.difmise1, aes(x = log10(n), y = log10(hh.sq))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
    # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$sig.sqdiff)), 
    #          label = paste0("Slope: ", round(slope_sig, 2)), hjust = 0, color = "red",vjust = 1) +
    theme_bw() + xlab("log10(n)") + ylab("log10(HH)")
  
  p7 <- ggplot(keep.difmise1, aes(x = log10(n), y = log10(hh.sqp))) +
    geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
    # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$sig.sqdiff)), 
    #          label = paste0("Slope: ", round(slope_sig, 2)), hjust = 0, color = "red",vjust = 1) +
    theme_bw() + xlab("log10(n)") + ylab("log10(hh)")
  
  # Arrange plots
  grid.arrange(p1, p2, p3, p4,p5,p6,p7, ncol = 4, nrow = 2)
  
  Sys.sleep(1)  # Adds a slight delay to allow visualization
}

###### Investigating the behavior of Atilde #####
## function for infinity norm
infinity_norm_matrix <- function(A) {
  max(rowSums(abs(A)))
}
library(ggplot2)
library(gridExtra)
library(dplyr)
r=1
# nseg = 37
pord = 1
bdeg = 3
sig <- 0.1
# nsim <- 10
lambda <-2
# Define the sample sizes and smoothing parameters
log.n <- seq(log10(10000), log10(100000), log10(1.11))  # Logarithmic grid for n
ns <- ceiling(10^log.n)  # Convert log scale to actual sample sizes
v <- (2 * (2 * pord)) + 1  # Smoothness parameter
rate <- 1 / v
Ks <- ceiling(10^(1.1 + rate * log.n))  # Number of segments

keep.df<- data.frame(
  n = rep(NA,length(ns)),
  K = rep(NA, length(ns)),
  L_inf.Atilde= rep(NA, length(ns)),
  L_2.Atilde.x= rep(NA, length(ns)),
  L_inf.AT= rep(NA, length(ns)),
  L_2.AT.x = rep(NA, length(ns))
)
# Define the specific value of x you want to include
target_x <- 0.5

for (j in 1:length(ns)) {

  # Generate the sequence for a given sample size ns[j]
  x <- sort((c(seq(0, 1, length.out = (ns[j]-1) ), target_x)))
 length(x)
  f <- 32 * exp(-8 * (1 - 2 * x)^2) * (1 - 2 * x)
  # if (r == 1) {
  #   fr <- (4096 * x^2 - 4096 * x + 960) * exp(-8 * (1 - 2 * x)^2)
  # } else if (r == 2) {
  #   fr <- -(262144 * x^3 - 393216 * x^2 + 184320 * x - 26624) * exp(-8 * (1 - 2 * x)^2)
  # }
    y <- f + rnorm(n = ns[j], mean = 0, sd = sig)
    length(y)
    beh <- AtildeA(x=x,y=y,lambda=lambda, r=r, x.grid=x, nseg=Ks[j]
                   ,pord=pord,bdeg=bdeg)
    idx <- which(x == target_x)
    Atilde <- beh$Atilde
    AT<- beh$AT
    dim(Atilde)
    dim(AT)
    
   Atilde.L_infty.norm <-  infinity_norm_matrix(Atilde)
   Atilde.x <-as.vector(Atilde[idx,])
   Atilde.L_2.norm <- max(Atilde.x)
     # norm(x=Atilde.x, type = "2")
   print(Atilde.L_2.norm)
  Sys.sleep(2)
   AT.L_infty.norm <-  infinity_norm_matrix(AT)
   AT.x <-as.vector(AT[floor(Ks[j]/2)],)
   AT.L_2.norm <- abs(max(AT.x))
     # norm(x=AT.x, type = "2")
   
    keep.df$n[j] <- ns[j]
    keep.df$K[j] <- Ks[j]
    keep.df$L_inf.Atilde[j] <- Atilde.L_infty.norm
    keep.df$L_2.Atilde.x[j]<- Atilde.L_2.norm
    keep.df$L_inf.AT[j] <- AT.L_infty.norm
    keep.df$L_2.AT.x[j]<- AT.L_2.norm
    # counter <- counter + 1
    print(j)
    print(ns[j])
    print(x[idx])
    keep.df1 <- keep.df %>%
      filter(!is.na(L_inf.Atilde), !is.na(L_2.Atilde.x), 
             !is.na(L_inf.AT), !is.na(L_2.AT.x))
    # Compute slopes
    slope_L_inf.Atilde <- coef(lm(log10(L_inf.Atilde) ~ log10(n), data = keep.df1))[2]
    slope_L_2.Atilde.x <- coef(lm(log10(L_2.Atilde.x) ~ log10(n), data = keep.df1))[2]
    slope_L_inf.AT <- coef(lm(log10(L_inf.AT) ~ log10(n), data = keep.df1))[2]
    slope_L_2.AT.x <- coef(lm(log10(L_2.AT.x) ~ log10(n), data = keep.df1))[2]
    
    # Print slopes in console for debugging
    print(paste("Slope of L_inf.Atilde:", round(slope_L_inf.Atilde, 3)))
    print(paste("Slope of L_2.Atilde.x:", round(slope_L_2.Atilde.x, 3)))
    print(paste("Slope of L_inf.AT:", round(slope_L_inf.AT, 3)))
    print(paste("Slope of L_2.AT.x:", round(slope_L_2.AT.x, 3)))
    
    # Plotting with slopes
    p1 <- ggplot(keep.df1, aes(x = log10(n), y = log10(L_inf.Atilde))) +
      geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue")+ 
      # + 
      #   annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$mise.diff)), 
      #            label = paste0("Slope: ", round(slope_mise, 2)), hjust = 0, color = "red",vjust = 1) +
      theme_bw() + xlab("log10(n)") + ylab("log10(L_inf.Atilde)")
    
    p2 <- ggplot(keep.df1, aes(x = log10(n), y = log10(L_2.Atilde.x))) +
      geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
      # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$var.diff)), 
      #          label = paste0("Slope: ", round(slope_var, 2)), hjust = 0, color = "red",vjust = 1) +
      theme_bw() + xlab("log10(n)") + ylab("log10(L_2.Atilde.x)")
    
    p3 <- ggplot(keep.df1, aes(x = log10(n), y = log10(L_inf.AT))) +
      geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
      # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$sq.bias.diff)), 
      #          label = paste0("Slope: ", round(slope_bias, 2)), hjust = 0,vjust = 1, color = "red", size = 5) +
      theme_bw() + xlab("log10(n)") + ylab("log10(L_inf.AT)")
    
    p4 <- ggplot(keep.df1, aes(x = log10(n), y = log10(L_2.AT.x))) +
      geom_point() + geom_smooth(method = "lm", se = FALSE, color = "blue") + 
      # annotate("text", x = min(log10(keep.difmise1$n)), y = max(log10(keep.difmise1$sig.sqdiff)), 
      #          label = paste0("Slope: ", round(slope_sig, 2)), hjust = 0, color = "red",vjust = 1) +
      theme_bw() + xlab("log10(n)") + ylab("log10(L_2.AT.x)")
    
    # Arrange plots
    grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
    
    Sys.sleep(1)  # Adds a slight delay to allow visualization
    
}

beh <- AtildeA(x=x,y=y,lambda=0.1, r=0, x.grid, nseg=35,pord=2,bdeg=3)