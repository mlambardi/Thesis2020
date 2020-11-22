require(dplyr)
require(data.table)

dir.create("~/Temp/tsn/data", recursive = T)

# require(sn) # skew-normal distribution

seed <- 1234
dofs <- c(10, 20, 30, 40, 50) # t-Student distro, hp: dof > 2
alphas <- c(0, 5, 10, 20) # skew-normal, values as in Sartori (2006)
ds <- c(10, 20, 50) # number of predictors
noverds <- c(10, 20, 50)
trials <- 10000

set.seed(seed)

maxd <- max(ds)
maxn <- maxd*max(noverds)

save(seed, dofs, alphas, ds, noverds, trials, file="configs.RData")

X <- cbind(1, matrix(rnorm((maxd-1)*maxn), nrow = maxn, ncol = maxd-1))
colnames(X) <- paste0("X",1:maxd)
X <- as.matrix(X)

configs <- expand.grid(
  d=ds,
  noverd=noverds,
  model=c(paste0("t", dofs), paste0("sn", alphas))
)
configs$theta <- as.numeric(sub("^[a-z]+(\\d.*)$", "\\1", configs$model))
configs$model <- as.character(sub("^([a-z]+)\\d.*$", "\\1", configs$model))
configs$n <- configs$d*configs$noverd
configs <- lapply(1:nrow(configs), function(i) {
  need <- configs$n[i]*trials/min(configs$n*trials)
  cbind(configs[rep(i, need),], index=1:need)
}) %>%
  data.table::rbindlist() %>%
  as.data.frame()
configs$trials <- trials*min(configs$n)/configs$n

pbapply::pbapply(configs, 1, function(v) {
  trials <- as.numeric(v["trials"])
  d <- as.numeric(v["d"])
  n <- as.numeric(v["n"])
  model <- as.character(v["model"])
  theta <- as.numeric(v["theta"])
  index <- as.numeric(v["index"])
  
  Z <- matrix(switch(model,
     t=rt(n*trials, df = theta),
     sn=sn::rsn(n*trials, alpha = theta)
  ), nrow = n)
  aux <- X[1:n,1:d,drop=F]
  Z <- t(Z - aux %*% (solve(crossprod(aux)) %*% crossprod(aux, Z)))
  Z <- Z / sqrt(rowSums(Z^2)/(n - d))
  save(n, d, theta, model, Z, file=paste0("~/Temp/tsn/data/", model, theta, "d", d, "n", n, "take", index, ".RData"))
})
