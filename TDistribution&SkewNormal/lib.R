
# theta is the shape parameter, so
# * d.o.f. in the t case
# * alpha in the skew-normal case, such that the implied skewness is
#   0.5*(4 - pi)*mu^3/(1 - mu^2)^1.5, with mu=sqrt(2/pi)*alpha/sqrt(1 + alpha^2)
# z is a vector of linear model residuals, standardized with respect to the square root of d.o.f.-corrected variance
# if z is a matrix, each row is treated as a distinct dataset, so that rows are independent

models <- list(
  general=list( # just a wrapper. model can be any element of models but the general one
    name="interface",
    description="general class that deparses models into usable objects",
    score=function(theta, z, d, beta, sigma, model) {
      model$qt(theta = theta, z = z) - models$general$scoreadj(theta = theta, z = z, d = d, beta = beta, sigma = sigma, model = model)
    },
    se=function(theta, z, model) {
      n <- ncol(z)
      qt <- model$qt(theta = theta, z = z)
      q_t_t <- rowMeans(qt^2)
      qtt <- rowMeans(model$qtt(theta = theta, z = z))
      qtz <- model$qtz(theta = theta, z = z)
      qtz0 <- rowMeans(qtz)
      qtz1 <- rowMeans(qtz*z)
      qt1 <- rowMeans(qt*z)
      qt2 <- rowMeans(qt*z^2)
      rho3 <- model$rho(k = 3, theta = theta)
      rho4 <- model$rho(k = 4, theta = theta)
      sqrt((qtz0^2 + qtz0*qtz1*rho3 - 2*qtz0*qt1 + qtz1^2*(rho4 - 1)/4 - qtz1*qt2 + q_t_t)/n)/abs(qtt)
    },
    scoreadj=function(theta, z, d, beta, sigma, model) {
      if (!beta & !sigma) {
        return(0)
      }
      n <- ncol(z)
      if (model$exact) {
        qtz1 <- model$qtz_k(theta = theta, k = 1) # it's needed anyways in NC estimation
      } else {
        qtz <- model$qtz(theta = theta, z = z)
        qtzz <- model$qtzz(theta = theta, z = z)
        qtz1 <- qtz*z
      }
      if (beta) {
        if (model$exact) {
          qtzz0 <- model$qtzz_k(theta = theta, k = 0)
        } else {
          qtzz0 <- qtzz
        }
        adjb <- d*(qtzz0/2 - qtz1)
      } else {
        adjb <- 0
      }
      if (sigma) {
        if (model$exact) {
          qtz3 <- model$qtz_k(theta = theta, k = 3)
          qtzz2 <- model$qtzz_k(theta = theta, k = 2)
        } else {
          qtz3 <- qtz*z^3
          qtzz2 <- qtzz*z^2
        }
        rho4 <- model$rho(k = 4, theta = theta)
        adjs <- (qtz1 - qtz3)/2 + (qtzz2 + 3*qtz1)*(rho4 - 1)/8
      } else {
        adjs <- 0
      }
      if (beta & sigma) {
        if (model$exact) {
          qtz0 <- model$qtz_k(theta = theta, k = 0)
          qtzz1 <- model$qtzz_k(theta = theta, k = 1)
        } else {
          qtz0 <- qtz
          qtzz1 <- qtzz*z
        }
        rho3 <- model$rho(k = 3, theta = theta)
        adjbs <- (qtz0 + qtzz1)*rho3/2
      } else {
        adjbs <- 0
      }
      return((adjb + adjbs + adjs)/n)
    }
  ),
  t=list(
    name="t",
    description="t distribution with theta degrees of freedom, standardized",
    exact=T,
    convert=function(eta, beta, sigma) { # 0 < eta < 1, function must be increasing
      ifelse(
        eta==0,
        2*(1 + sigma),
        ifelse(
          eta==1,
          +Inf,
          2*(1 + sigma)/abs(1 - eta) # 2 or 4 < result < +Inf
        )
      )
    },
    rho=function(k, theta) {
      ifelse(
        (theta > k) & (k >= 0),
        if (k %% 2 == 0) {
          if (k <= 2) {
            1
          } else {
            aux <- 1
            for (i in 1:round(k/2)) {
              aux <- aux*(2*i - 1)*(theta - 2)/(theta - 2*i)
            }
            aux
          }
        } else {
          0
        },
        NA
      )
    },
    llik=function(theta, z, w=1/(1 + z^2/(theta - 2))) {
      log(w)*(theta + 1)/2 + models$t$lconst(theta = theta, deriv = 0)
    },
    qt=function(theta, z, w=1/(1 + z^2/(theta - 2))) { # score function
      (log(w) + w*z^2*(theta + 1)/(theta - 2)^2)/2 + models$t$lconst(theta = theta, deriv = 1)
    },
    qtt=function(theta, z, w=1/(1 + z^2/(theta - 2))) {
      (w^2*z^4*(theta + 1)/(theta - 2)^4 - 6*w*z^2/(theta - 2)^3)/2 + models$t$lconst(theta = theta, deriv = 2)
    },
    qtz=function(theta, z, w=1/(1 + z^2/(theta - 2))) {
      3*w*z/(theta - 2)^2 - w^2*z^3*(theta + 1)/(theta - 2)^3
    },
    qtz_k=function(theta, k) {
      3*models$t$wz(theta=theta, w=1, z=1+k)/(theta - 2)^2 - models$t$wz(theta=theta, w=2, z=3+k)*(theta + 1)/(theta - 2)^3
    },
    qtzz=function(theta, z, w=1/(1 + z^2/(theta - 2))) {
      3*w/(theta - 2)^2 - 3*w^2*z^2*(theta + 3)/(theta - 2)^3 + 4*w^3*z^4*(theta + 1)/(theta - 2)^4
    },
    qtzz_k=function(theta, k) {
      3*models$t$wz(theta=theta, w=1, z=k)/(theta - 2)^2 - 3*models$t$wz(theta=theta, w=2, z=2+k)*(theta + 3)/(theta - 2)^3 + 4*models$t$wz(theta=theta, w=3, z=4+k)*(theta + 1)/(theta - 2)^4
    },
    wz=function(theta, w, z) { # E[W^w * Z^z;theta]
      models$t$rho(theta=theta + 2*w, k=z)*exp(
        models$t$lconst(theta=theta, deriv=0) - models$t$lconst(theta=theta + 2*w, deriv=0)
      )*((theta - 2)/(theta + 2*w - 2))^((z + 1)/2)
    },
    lconst=function(theta, deriv=0) {
      mylgamma <- function(x, deriv) {
        if (deriv==0) {
          lgamma(x)
        } else {
          psigamma(x = x, deriv = deriv - 1)
        }
      }
      auxlog <- function(x, deriv) {
        (if (deriv==0) {
          -log(pi*x)
        } else {
          (-1/x)^deriv * factorial(deriv-1)
        })/2
      }
      # returns the deriv-th derivative of the log-normalization constant
      ifelse(
        theta <= 2,
        NA,
        (mylgamma(x = (theta + 1)/2, deriv = deriv) - mylgamma(x = theta/2, deriv = deriv))/2^deriv + auxlog(x = theta - 2, deriv = deriv)
      )
    }
  ),
  sn=list(
    name="sn",
    description="skew-normal distribution with skewness parameter theta, standardized",
    exact=F,
    convert=function(eta, beta, sigma) {
      ifelse(
        eta==0,
        -Inf,
        ifelse(
          eta==1,
          +Inf,
          tan(pi*(eta - 1/2))
        )
      )
    },
    rho=function(k, theta, b=sqrt(2/pi)*theta/sqrt(1 + theta^2)) {
      if (k==3) {
        (4 - pi)/2*b^3/(1 - b^2)^(3/2)
      } else {
        if (k==4) {
          3 + 2*(pi - 3)*b^4/(1 - b^2)^2
        } else {
          stop("unimplemented, not needed anyways")
        }
      }
    },
    provide=function(theta, Z) {
      # All uppercase quantities are matrices in general
      # All lowercase quantities are vectors or constant
      # Matrices must be first in products and sums
      b <- sqrt(2/pi)*theta/sqrt(1 + theta^2)
      a <- sqrt(1 - b^2)
      a_b <- -b/a
      b_t <- sqrt(2/pi)*(1 + theta^2)^(-3/2)
      a_t <- a_b*b_t
      b_tt <- -3*sqrt(2/pi)*theta*(1 + theta^2)^(-5/2)
      a_bb <- -(a_b^2 + 1)/a
      a_tt <- a_bb*b_t^2 + a_b*b_tt
      w_z <- a # beware, must be last in products, since vector
      w_zz <- 0
      w_tz <- a_t # as well
      w_tzz <- 0 # as well
      W <- Z*a + b
      W_t <- Z*a_t + b_t
      W_tt <- Z*a_tt + b_tt
      return(list(a=a, a_t=a_t, a_tt=a_tt, W=W, W_t=W_t, W_tt=W_tt, w_z=w_z, w_zz=w_zz, w_tz=w_tz, w_tzz=w_tzz))
    },
    llik=function(theta, z, dat=models$sn$provide(theta = theta, Z = z)) {
      with(dat, -W^2/2 + pnorm(W*theta, log.p = T) + log(a)) # why log.p instead of log?
    },
    qt=function(theta, z, dat=models$sn$provide(theta = theta, Z = z)) {
      with(
        dat,
        {
          M0 <- models$sn$millsratio(x = W*theta, deriv = 0)
          -W*W_t + M0*(W + W_t*theta) + a_t/a
        }
      )
    },
    qtt=function(theta, z, dat=models$sn$provide(theta = theta, Z = z)) {
      with(
        dat,
        {
          M0 <- models$sn$millsratio(x = W*theta, deriv = 0)
          M1 <- models$sn$millsratio(x = W*theta, deriv = 1)
          -(W_t^2 + W*W_tt) + M1*(W + W_t*theta)^2 + M0*(W_t*2 + W_tt*theta) + (a_tt*a - a_t^2)/a
        }
      )
    },
    qtz=function(theta, z, dat=models$sn$provide(theta = theta, Z = z)) {
      with(
        dat,
        {
          M0 <- models$sn$millsratio(x = W*theta, deriv = 0)
          M1 <- models$sn$millsratio(x = W*theta, deriv = 1)
          -(W_t*w_z + W*w_tz) + M1*(W + W_t*theta)*w_z*theta + M0*(w_z + w_tz*theta)
        }
      )
    },
    qtzz=function(theta, z, dat=models$sn$provide(theta = theta, Z = z)) {
      with(
        dat,
        {
          M0 <- models$sn$millsratio(x = W*theta, deriv = 0)
          M1 <- models$sn$millsratio(x = W*theta, deriv = 1)
          M2 <- models$sn$millsratio(x = W*theta, deriv = 2)
          -(W_t*w_zz + W*w_tzz + 2*w_z*w_tz) + ((W + W_t*theta)*(M2*theta*w_z^2 + M1*w_zz)*theta + M0*(w_zz + w_tzz*theta) + 2*M1*theta*w_z*(w_z + w_tz*theta))
        }
      )
    },
    millsratio=function(x, deriv=0) {
      mr0 <- function (x) { # after VGAM::mills.ratio, but error-safe
        ifelse(x < -100, -x/(1 - 1/x^2 + 3/x^4), exp(dnorm(x, log=T) - pnorm(x, log.p=T)))
      }
      mr1 <- function(x, m0=mr0(x)) {
        -m0*(x + m0)
      }
      mr2 <- function(x, m0=mr0(x), m1=mr1(x)) {
        m0*((x + m0)*(x + 2*m0) - 1)
      }
      if (deriv == 0) {
        mr0(x)
      } else {
        if (deriv == 1) {
          mr1(x)
        } else {
          if (deriv == 2) {
            mr2(x)
          } else {
            stop("unimplemented, though not needed")
          }
        }
      }
    }
  ),
  gn=list(
    name="gn",
    description="exponential power/generalized error/generalized normal",
    exact=F,
    convert=function(eta, beta, sigma) {
      ifelse(eta==0, 2, ifelse(eta==1, +Inf, 2/abs(1 - eta)))
    },
    rho=function(k, theta, abs=F) { # 
      if (abs) {
        exp(lgamma((k + 1)/theta) + lgamma(1/theta)*(k/2 - 1) - lgamma(3/theta)*k/2)
      } else {
        ifelse(k %% 2 == 1, 0, models$gn$rho(k = k, theta = theta, abs = T))
      }
    },
    llik=function(theta, z) {
      ls <- models$gn$logstdev(theta = theta, deriv = 0)
      lc <- models$gn$lconst(theta = theta, deriv = 0)
      -abs(z)^theta*exp(-theta*ls) + lc - ls
    },
    qt=function(theta, z) {
      ls0 <- models$gn$logstdev(theta = theta, deriv = 0)
      ls1 <- models$gn$logstdev(theta = theta, deriv = 1)
      lc1 <- models$gn$lconst(theta = theta, deriv = 1)
      -abs(z)^theta*exp(-theta*ls0)*(log(z) - ls0 - theta*ls1) + lc1 - ls1
    },
    qtt=NULL,
    qtz=NULL,
    qtz_k=NULL,
    qtzz=NULL,
    qtzz_k=NULL,
    logstdev=function(theta, deriv=0) {
      if (deriv == 0) {
        lgamma(3/theta) - lgamma(1/theta)
      } else {
        if (deriv == 1) {
          -(digamma(3/theta)*3 - digamma(1/theta))/theta^2
        } else { # deriv == 2
          (trigamma(3/theta)*9 - trigamma(1/theta))/theta^4 + 2*(digamma(3/theta)*3 - digamma(1/theta))/theta^3
        }
      }/2
    },
    lconst=function(theta, deriv=0) {
      if (deriv == 0) {
        log(theta/2) - lgamma(1/theta)
      } else {
        if (deriv == 1) {
          1/theta + digamma(1/theta)/theta^2
        } else { # deriv == 2
          -1/theta^2 - 2*digamma(1/theta)/theta^3 - trigamma(1/theta)/theta^4
        }
      }
    }
  )
)

# estimate a single shape parameter under linear models based on residuals
# residuals must be standardized with respect to est sigma corrected for d.o.f
estResModel <- function(Z, d, mod, corrbeta=F, corrsigma=F, enhancesigma=F, verbose=F, eps=1e-5) {
  minest <- rep(0, nrow(Z))
  est <- rep(NA, nrow(Z))
  maxest <- rep(1, nrow(Z))
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    print(setTxtProgressBar(pb, 0))
  }
  index <- 1:nrow(Z)
  # aux <- matrix(NA, ncol = 0, nrow = length(finite))
  iters <- ceiling(log(eps)/log(0.5)) + 10 # makes it possible to actually reach the boundary estimates
  if (!enhancesigma) {
    Z <- Z/sqrt(1 - d/ncol(Z))
  }
  for (i in 1:iters) {
    est[index] <- (maxest[index] + minest[index])/2
    aux1 <- (est[index] < eps)
    aux2 <- (est[index] > 1-eps)
    est[index[aux1]] <- 0
    est[index[aux2]] <- 1
    index <- index[!aux1 & !aux2]
    if (length(index) == 0) break
    # aux <- cbind(aux, est[finite])
    pt <- mod$convert(eta = est[index], beta = corrbeta, sigma = corrsigma)
    score <- rowSums(models$general$score(theta = pt, z = Z[index,,drop=F], d = d, beta = corrbeta, sigma = corrsigma, model = mod))
    minest[index] <- ifelse(score >= 0, est[index], minest[index])
    maxest[index] <- ifelse(score <= 0, est[index], maxest[index])
    if (verbose) {
      print(setTxtProgressBar(pb, i/iters))
    }
  }
  if (verbose) {
    print(setTxtProgressBar(pb, 1))
    close(pb)
  }
  mod$convert(eta = est, beta = corrbeta, sigma = corrsigma)
}

# models$t$exact <- T
# 
# res <- data.frame(
#   bs=estResModel(Z = Z, d = d, mod = models$t, corrbeta = F, corrsigma = F, enhancesigma = F, verbose = F),
#   bS=estResModel(Z = Z, d = d, mod = models$t, corrbeta = F, corrsigma = T, enhancesigma = T, verbose = F),
#   Bs=estResModel(Z = Z, d = d, mod = models$t, corrbeta = T, corrsigma = F, enhancesigma = T, verbose = F),
#   BS=estResModel(Z = Z, d = d, mod = models$t, corrbeta = T, corrsigma = T, enhancesigma = T, verbose = F)
# )
# 
# require(ggplot2)
# require(dplyr)
# 
# type <- geom_violin
# 
# (log(log(res-4+1)) - log(log(theta-4+1))) %>%
#   ggplot() +
#   type(aes(x=bs, y="bs")) +
#   type(aes(x=bS, y="bS")) +
#   type(aes(x=Bs, y="Bs")) +
#   type(aes(x=BS, y="BS"))
