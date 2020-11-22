
dir.create("~/Temp/tsn/rho", recursive = T)

# Z is a vector of regression residuals, standardized with respect to the d.o.f.-corrected estimate of standard deviation
# if Z is a matrix, each row is a distinct dataset, all datasets being iid replicates of each other
# d is the number of regressors, the model is assumed to include the intercept, so that c/n=1
# M is the maximum order of standardized moments to be estimated
# maxiter is the maximum number of iterations in three-blocks nuisance control estimation
# tol is the maximum absolute value of relative update in the three-blocks nuisance control estimation to assert convergence
standardized.moments <- function(Z, d, M=6, maxiter=30, tol=1e-5) {
  if (is.vector(Z)) {
    Z <- matrix(Z, nrow = 1)
  }
  M <- max(M, 6)
  n <- ncol(Z)
  trials <- nrow(Z)
  
  rho <- list()
  
  ks <- 3:M
  rho$naive <- matrix(sapply(ks, function(k) rowMeans(Z^k)), ncol = length(ks))
  colnames(rho$naive) <- paste0("rho", ks)
  
  # P estimation, standardizing with respect to biased variance estimator
  rho$pe <- rho$naive %*% diag(1/(1 - d/n)^(ks/2))
  colnames(rho$pe) <- colnames(rho$naive)
  
  # NC estimation
  rho$nc2 <- rho$pe
  rho$nc2[,"rho3"] <- rho$naive[,"rho3"]/(1 - 3*d/n)
  rho$nc2[,"rho4"] <- (rho$naive[,"rho4"] - 6*d/n)/(1 - 4*d/n)
  for (k in 5:M) {
    rho$nc2[,paste0("rho", k)] <- if (n > k*d) {
      rk_naive <- rho$naive[,paste0("rho", k)]
      rkm2_nc2 <- rho$nc2[,paste0("rho", k-2)]
      (rk_naive - choose(k, 2)*d/n*rkm2_nc2)/(1 - k*d/n)
    } else {
      warning(paste0("insufficient data, ", k, "-th order moments will not be estimated under two-blocks nuisance control"))
      NA
    }
  }
  rho$nc3 <- cbind(rho1=0, rho2=1, rho$naive)
  k1 <- 3:(M - 2)
  k2 <- (M - 1):M
  corrbeta <- T
  corrsigma <- T
  for (i in 1:maxiter) {
    rkm2 <- rho$nc3[,k1-2,drop=F]
    rkm1 <- rho$nc3[,k1-1,drop=F]
    rk <- rho$nc3[,k1,drop=F]
    rkp2 <- rho$nc3[,k1+2,drop=F]
    r3 <- rho$nc3[,3]
    r4 <- rho$nc3[,4]
    Ab <- d*(rkm2 %*% diag(choose(k1, 2)) - rk %*% diag(k1))
    Ab2 <- d*(rho$nc3[,k2-2,drop=F] %*% diag(choose(k2, 2)) - rho$nc3[,k2,drop=F] %*% diag(k2))
    As <- (rk - rkp2) %*% diag(k1) / 2 + rk * tcrossprod(r4 - 1, k1^2 + 2*k1) / 8
    Abs <- rkm1 * tcrossprod(r3, k1^2) / 2
    upd <- rho$naive - rho$nc3[,3:M] - cbind(
      Ab + Abs + As,
      Ab2
    )/n
    rho$nc3[,3:M] <- rho$nc3[,3:M] + upd
    # rho$nc3[,4] <- pmax(rho$nc3[,4], 1) # this is being forced
    if (i == 1) {
      rho$nc3onestep <- rho$nc3[,3:M,drop=F]
    }
  }
  rho$nc3 <- rho$nc3[,-c(1:2),drop=F]
  colnames(upd) <- colnames(rho$nc3onestep) <- colnames(rho$nc3) <- colnames(rho$pe)
  rho$nc3 <- ifelse(abs(upd/rho$nc3) < tol, rho$nc3, NA)
  
  for (m in names(rho)) {
    colnames(rho[[m]]) <- paste0(colnames(rho[[m]]), m)
  }
  rho <- do.call("cbind", rho)
  colnames(rho) <- sub("^.*?rho", "rho", colnames(rho))
  return(rho)
}

Zes <- t(apply(matrix(rt(100*1000, 20), nrow = 1000), 1, scale, center = T, scale = T))
summary(standardized.moments(
  Z=Zes,
  d=1, M=5
)[,c("rho4pe", "rho4naive", "rho4nc2", "rho4nc3")]-(3 + 6/(20-4)))

pbapply::pblapply(list.files(path = "~/Temp/tsn/data", full.names = T), function(f) {
  load(f)
  iters <- 30
  tol <- 1e-5
  M <- 8
  tryCatch({
    rho <- standardized.moments(Z, d=d, M=M, maxiter=iters, tol=tol)
  }, error=function(e) {
    save(Z, d, M, iters, file="dump.RData")
    stop("dumped")
  })
  save(model, theta, d, n, M, rho, iters, tol, file=sub("/Temp/tsn/data/", "/Temp/tsn/rho/", f))
})
