######## Comparing with other methods ##########################################

# Load required libraries
library(sfsmisc)
library(locpol)

####### Charnigo functions #######################
### Functions for use in multiple simulation studies

# Function to compute the C matrix
Cmatrix <- function(x, weights) {
  n <- length(x)
  result <- matrix(0, n, n)
  result[1, ] <- c(-1 / (x[2] - x[1]), 1 / (x[2] - x[1]), rep(0, n - 2))
  result[n, ] <- c(rep(0, n - 2), -1 / (x[n] - x[n - 1]), 1 / (x[n] - x[n - 1]))

  for (i in 2:(n - 1)) {
    j <- min(i, n - i + 1, length(weights) + 1) - 1
    gaps <- rep(0, j)
    for (k in 1:j) {
      gaps[k] <- x[i + k] - x[i - k]
    }
    weightsn <- weights[1:j] / sum(weights[1:j]) / gaps
    result[i, ] <- c(rep(0, i - j - 1), -weightsn[j:1], 0, weightsn[1:j], rep(0, n - i - j))
  }
  return(result)
}

# Function to compute the C2 matrix
C2matrix <- function(x, weights1, weights2) {
  A <- Cmatrix(x, weights1)
  B <- Cmatrix(x, weights2)
  C <- B %*% A
  return(C)
}

# Function for first derivative estimation
nedq1no <- function(i, y, x, weights) {
  n <- length(x)
  result <- 0
  if (i == 1) {
    result <- (y[2] - y[1]) / (x[2] - x[1])
  } else if (i == n) {
    result <- (y[n] - y[n - 1]) / (x[n] - x[n - 1])
  } else {
    j <- min(i, n - i + 1, length(weights) + 1) - 1
    weightsn <- weights[1:j] / sum(weights[1:j])
    for (k in 1:length(weightsn)) {
      result <- result + weightsn[k] * (y[i + k] - y[i - k]) / (x[i + k] - x[i - k])
    }
  }
  return(result)
}

# Function for second derivative estimation
nedq2no <- function(i, y, x, weights) {
  n <- length(x)
  result <- 0
  if (i < 3) {
    result <- (-(y[2] - y[1]) / (x[2] - x[1]) + (y[3] - y[2]) / (x[3] - x[2])) / (x[2] - x[1])
  } else if (i > n - 2) {
    result <- (-(y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) + (y[n] - y[n - 1]) / (x[n] - x[n - 1])) / (x[n] - x[n - 1])
  } else {
    j <- min(floor((i + 1) / 2), floor((n - i + 2) / 2), length(weights) + 1) - 1
    weightsn <- weights[1:j] / sum(weights[1:j])
    for (k in 1:length(weightsn)) {
      result <- result + weightsn[k] * ((y[i + 2 * k] - y[i]) / (x[i + 2 * k] - x[i]) - (y[i] - y[i - 2 * k]) / (x[i] - x[i - 2 * k])) / (x[i + k] - x[i - k])
    }
  }
  return(result)
}

# Functions to generate weights
weights1 <- function(k) {
  w1 <- rep(0, k)
  sum1 <- 0
  for (i in 1:k) {
    sum1 <- sum1 + i^2
  }
  for (i in 1:k) {
    w1[i] <- i^2 / sum1
  }
  return(w1)
}

weights2 <- function(k) {
  w2 <- rep(0, k)
  sum1 <- 0
  for (i in 1:k) {
    sum1 <- sum1 + i
  }
  for (i in 1:k) {
    w2[i] <- i / sum1
  }
  return(w2)
}

# Function to find the minimum argument in a matrix
matrixargmin <- function(matrix) {
  norow <- dim(matrix)[1]
  nocol <- dim(matrix)[2]
  temp <- which.min(matrix)
  colargmin <- as.integer((temp + norow - 1) / norow)
  rowargmin <- temp - (colargmin - 1) * norow
  return(rowargmin, colargmin)
}

# Function to compute L matrix
Lmatrix <- function(x, h, nknots = 35) {
  n <- length(x)
  Lmatrix <- matrix(0, n, n)
  for (k in 1:n) {
    tempy <- rep(0, n)
    tempy[k] <- 1
    Lmatrix[, k] <- predict(smooth.spline(x = x, y = tempy, spar = h, nknots = nknots), x, deriv = 0)$y
  }
  lambda <- smooth.spline(x = x, y = tempy, spar = h, nknots = 9)$lambda
  return(list(Lmatrix = Lmatrix, lambda = lambda))
}

# Function to compute L1 matrix
L1matrix <- function(x = x, h = -0.3, x.grid = x.grid, nknots = 35) {
  n <- length(x)
  L1matrix <- matrix(0, n, n)
  for (k in 1:n) {
    tempy <- rep(0, n)
    tempy[k] <- 1
    L1matrix[, k] <- predict(smooth.spline(x = x, y = tempy, spar = h, nknots = nknots), x.grid, deriv = 1)$y
  }
  return(L1matrix)
}

# Function to compute L2 matrix
L2matrix <- function(x, h, x.grid = x.grid, nknots = 35) {
  n <- length(x)
  L2matrix <- matrix(0, n, n)
  for (k in 1:n) {
    tempy <- rep(0, n)
    tempy[k] <- 1
    L2matrix[, k] <- predict(smooth.spline(x = x, y = tempy, spar = h, nknots = nknots), x.grid, deriv = 2)$y
  }
  return(L2matrix)
}

# Function to calculate derivative first estimate via GCp
GCP.est <- function(h = -0.23, x = x, y = y, k, nknots = 35) {
  n <- length(y)
  k1 <- k
  S1 <- c(rep(0, n * 0.05), rep(1, n - 2 * n * 0.05), rep(0, n * 0.05))

  L0Spline <- array(0, c(n, n, 1))
  L1Spline <- array(0, c(n, n, 1))
  L0Spline[, , 1] <- Lmatrix(x, h, nknots = nknots)$Lmatrix
  L1Spline[, , 1] <- L1matrix(x, h, x.grid = x, nknots = nknots)

  sigma2c <- sum((L0Spline[, , 1] %*% y - y)^2) / (n - sum(diag(L0Spline[, , 1])))

  hatmu1real <- L1Spline[, , 1] %*% y

  GCP1real <- array(0, length(k1))
  for (i in 1:length(k1)) {
    Emp1matrix <- array(0, c(n, n, length(k1[i]) + 1))
    Emp1matrix[, , 1] <- Cmatrix(x, weights1(k1[i]))
    Emp1matrix[, , 2] <- Cmatrix(x, weights1(1))

    empirical <- apply(as.matrix(c(1:length(x))), 1, nedq1no, y = y, x = x, weights = weights1(k1[i]))

    plot(x, empirical, type = "l")
    lines(x, hatmu1real, col = "red")

    GCP1 <- (c(hatmu1real) - empirical)^2 %*% S1 +
      sigma2c * sum(S1 %*% (2 * L1Spline[, , 1] * Emp1matrix[, , 1] - Emp1matrix[, , 1] * Emp1matrix[, , 1]))
    GCP1real[i] <- GCP1
  }

  loss <- mean(abs(GCP1real))

  return(loss)
}

# Function to calculate second derivative estimate via GCp
GCP.est2 <- function(h = 0.9, x = x, y = y, k1, k2, nknots = 35) {
  n <- length(y)
  k1 <- k1
  k2 <- k2
  S1 <- c(rep(0, n * 0.1), rep(1, n - 2 * n * 0.1), rep(0, n * 0.1))

  L0Spline <- array(0, c(n, n, 1))
  L2Spline <- array(0, c(n, n, 1))
  L0Spline[, , 1] <- Lmatrix(x, h, nknots = nknots)$Lmatrix
  L2Spline[, , 1] <- L2matrix(x, h, x.grid = x, nknots = nknots)

  sigma2c <- sum((L0Spline[, , 1] %*% y - y)^2) / (n - sum(diag(L0Spline[, , 1])))

  hatmu2real <- L2Spline[, , 1] %*% y

  plot(x, hatmu2real, type = "l")

  GCP2real <- array(0, length(k1) * length(k2))
  z <- 1
  for (i in 1:length(k1)) {
    for (j in 1:length(k2)) {
      Emp2matrix <- array(0, c(n, n, 2))
      Emp2matrix[, , 1] <- C2matrix(x, weights1(k1[i]), weights2(k2[j]))
      Emp2matrix[, , 2] <- C2matrix(x, weights1(1), weights2(1))

      empirical <- Emp2matrix[, , 1] %*% y
      GCP2 <- (c(hatmu2real) - c(empirical))^2 %*% S1 +
        sigma2c * sum(S1 %*% (2 * L2Spline[, , 1] * Emp2matrix[, , 1] - Emp2matrix[, , 1] * Emp2matrix[, , 1]))
      GCP2real[z] <- GCP2
      z <- z + 1
      plot(x, empirical, type = "l")
      lines(x, hatmu2real, col = "red")
    }
  }

  loss <- mean(abs(GCP2real))

  return(loss)
}

minh.fun2 <- function(hvec, x, y, x.grid, k1 = c(5, 15, 25), k2 = c(5, 15, 25), nknots = 35) {
  best <- sapply(hvec, function(h) {
    GCP.est2(h = h, x = x, y = y, k1 = k1, k2 = k2, nknots = nknots)
  })
  plot(hvec, best, type = "l")

  min_index <- which.min(best)
  opt.h <- hvec[min_index]
  L2Spline.2 <- L2matrix(x = x.grid, h = opt.h, x.grid = x.grid, nknots = nknots)
  fr.hat <- L2Spline.2 %*% y

  return(list(h = opt.h, fr.hat = fr.hat))
}

minh.fun <- function(hvec, x, y, x.grid = x.grid, k = c(5, 10, 15, 20, 25), nknots = 35) {
  best <- sapply(hvec, function(h) {
    GCP.est(h = h, x = x, y = y, k, nknots = nknots)
  })

  plot(hvec, best, type = "l")

  min_index <- which.min(best)
  opt.h <- hvec[min_index]
  L1Spline.1 <- L1matrix(x = x.grid, h = opt.h, x.grid = x.grid)
  fr.hat <- L1Spline.1 %*% y

  return(list(h = opt.h, fr.hat = fr.hat))
}

charnigo.est <- function(x = x, y = y, r = 1, x.grid, k = c(5, 10, 15, 20, 25), k1 = c(5, 15, 25),
                         k2 = c(5, 15, 25), nknots = 35, hs) {
  if (r == 1) {
    est.char <- minh.fun(hvec = hs, x, y, x.grid = x.grid, k = k, nknots = nknots)
    fr.hat <- est.char$fr.hat
  } else {
    est.char <- minh.fun2(hvec = hs, x, y, x.grid = x.grid, k1 = k1, k2 = k2, nknots = nknots)
    fr.hat <- est.char$fr.hat
  }

  return(list(x.grid = x.grid, h = est.char$h, fr.hat = fr.hat))
}

## End of Charnigo functions #######################





####### Dai et al. 2016 methods #######################

# Function to select smoothing parameters
bandwidth <- function(x, y, r0, p, q) {
  n <- length(x)
  m <- max(2, floor(n^(1 / 3)))
  sigma <- ord_r2(x, y, m)
  Int <- est_int_lop(x, y, q)

  MSE <- ord <- seq(floor(r0 / 2) - 3) * 0
  BiasC <- seq(floor(r0 / 2) - 3) * 0
  VarC <- seq(floor(r0 / 2) - 3) * 0
  k <- 1
  for (i in 4:floor(r0 / 2)) {
    r <- 2 * i
    I <- seq(2 * q + 1)
    I[1] <- r + 1
    for (ii in 2:(2 * q + 1)) {
      I[ii] <- sum((c(0:r) - r / 2)^(ii - 1))
    }

    temp <- Vmatrix(r, r / 2, q - 1)
    Var <- n^(2 * p) * temp[p + 1, p + 1] * abs(sigma)
    C1 <- 0
    for (j in 1:q) {
      C1 <- C1 + temp[j, p + 1] * I[q + j]
    }
    Bias <- (gamma(p + 1) / (gamma(q + 1) * n^(q - p)))^2 * C1^2 * abs(Int)
    BiasC[k] <- Bias
    VarC[k] <- Var
    MSE[k] <- Var + Bias
    ord[k] <- 2 * i
    k <- k + 1
  }
  seq_ord <- ord[order(MSE)[1]]
  bias_ord <- q

  index <- which.min(MSE)
  Bias_q <- BiasC[index]
  Var_q <- VarC[index]
  MSE_q <- min(MSE)

  return(list(params = c(MSE_q, bias_ord, seq_ord, MSE), Bias = Bias_q, Var = Var_q))
}

# Function to calculate d_pq
d_pq <- function(p, q, r, l, n) {
  d <- seq(r + 1)
  I <- seq(2 * q + 1)
  U <- matrix(0, q + 1, q + 1)
  I[1] <- r + 1
  for (i in 2:(2 * q + 1)) {
    I[i] <- sum((c(0:r) - l)^(i - 1))
  }
  for (i in 1:(q + 1)) {
    for (j in 1:(q + 1)) {
      U[i, j] <- I[i + j - 1]
    }
  }
  L <- seq(q + 1) * 0
  L[p + 1] <- -2 * gamma(p + 1) * n^p

  a <- solve(U, tol = 1e-300) %*% L
  for (k in 1:(1 + r)) {
    temp <- seq(q + 1)
    for (i in 1:(q + 1)) {
      temp[i] <- (k - 1 - l)^(i - 1)
    }
    d[k] <- -sum(a * temp) / 2
  }
  return(d)
}

# Function to estimate derivatives
der_opt <- function(x, y, order, q, r) {
  n <- length(x)
  est <- seq(n - r)
  d <- d_pq(order, q, r, floor(r / 2), n)
  for (i in (1 + floor(r / 2)):(n - (r - floor(r / 2)))) {
    est[i - floor(r / 2)] <- sum(d * y[(i - floor(r / 2)):(i + r - floor(r / 2))])
  }
  return(est)
}

# Function to estimate integral
est_int_lop <- function(x, y, q) {
  n <- length(x)
  d <- data.frame(x = x)
  d$y <- y
  r <- locpol(y ~ x, d, deg = q, xeval = d$x)
  int <- mean((r$lpFit[(0.1 * n + 1):(0.9 * n), q + 2])^2)
  return(int)
}

# Function to estimate optimal derivatives
est_opt <- function(x, y, order, q, r) {
  n <- length(y)
  D <- (max(x) - min(x)) * n / (n - 1)
  der_opt <- seq(n)
  d2 <- d_pq(order, q, r, r / 2, n) / D^order

  for (i in 1:n) {
    if (i <= (r / 2)) {
      d1 <- d_pq(order, q, r, i - 1, n) / D^order
      der_opt[i] <- sum(y[1:(r + 1)] * d1)
    }
    if (i >= (n - r / 2 + 1)) {
      d1 <- d_pq(order, q, r, i - n + r, n) / D^order
      der_opt[i] <- sum(y[(n - r):n] * d1)
    }
    if (i >= (r / 2 + 1) && i <= (n - r / 2)) {
      der_opt[i] <- sum(y[(i - r / 2):(i + r / 2)] * d2)
    }
  }

  return(der_opt)
}

# Function to calculate ord_r2
ord_r2 <- function(x, y, m) {
  n <- length(x)
  M <- n * m - m * (m + 1)
  dw <- ds <- fm <- fz <- 0
  s <- seq(length = m, from = 0, by = 0)
  d <- seq(length = m, from = 0, by = 0)
  w <- seq(length = m, from = 0, by = 0)
  for (j in 1:m) {
    w[j] <- (n - 2 * j) / M
    z1 <- y[1:(n - 2 * j)]
    z2 <- y[(1 + j):(n - j)]
    z3 <- y[(1 + 2 * j):(n)]
    z <- (z1 - 2 * z2 + z3)^2 / 6
    s[j] <- sum(z) / (n - 2 * j)
    d[j] <- (j / n)^2
  }
  dw <- sum(w * d)
  ds <- sum(w * s)
  fm <- sum(w * (d - dw)^2)
  fz <- sum(w * s * (d - dw))
  if (fm < 0.0000000000000001) {
    u <- ds
  } else {
    u <- ds - dw * fz / fm
  }
  return(u)
}

# Function to calculate V matrix
Vmatrix <- function(r, l, q) {
  I <- seq(2 * q + 1)
  U <- matrix(0, q + 1, q + 1)
  I[1] <- r + 1
  for (ii in 2:(2 * q + 1)) {
    I[ii] <- sum((c(0:r) - r / 2)^(ii - 1))
  }
  for (i in 1:(q + 1)) {
    for (j in 1:(q + 1)) {
      U[i, j] <- I[i + j - 1]
    }
  }

  a <- solve(U, tol = 1e-30)
  return(a)
}


dai.fit<- function(x, y,p=1,q=7,x.grid, fr.grid ){
  n<- length(x)
  r0=floor(n/2)
  res=bandwidth(x,y,r0=r0,p=p,q=q)
  first_der_est=est_opt(x.grid,y,order=p,q=res$params[2]-1,res$params[3])
  plot(x.grid,fr.grid,type="l")
  lines(x.grid,first_der_est,col='blue',pch=20)
  emise<- mean((first_der_est-fr.grid)^2)
  
  return(list(emise=emise, fr.hat=first_der_est ))
}
 ##### End of Dai et al. 2016 methods #######################


################## End of Functions for other methods ###########################




