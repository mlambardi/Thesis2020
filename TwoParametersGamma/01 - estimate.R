
# file.remove("~/Temp/gamma", recursive = T)
dir.create("~/Temp/gamma/data", recursive = T)

set.seed(1234)

alphas <- c(0.1,0.2,0.5,1,2,5,10)
l <- 1
ns <- c(5,10,20,50,100)
trials <- 100000
iters <- 30

source("../../codes/gamma/gammanew.R")

(pars <- c(n=min(ns),a=min(alphas),l=1))
sam <- matrix(rgamma(pars["n"]*10000, pars["a"], pars["l"]), ncol = pars["n"])
summary(gammaest(dat=sam, exact = F)[,c("ape.ye", "ape", "anc", "abr", "scorecritnc")])

pbapply::pbapply(as.matrix(expand.grid(n=ns, alpha=alphas)), 1, function(v) {
  a <- v[2]
  n <- v[1]
  dat <- matrix(rgamma(n*trials, a, l), nrow = trials)
  res <- gammaest(dat, iters = iters)
  save(a, l, n, trials, iters, res, file=paste0("~/Temp/gamma/data/gammaa", a, "n", n, ".RData"))
})
