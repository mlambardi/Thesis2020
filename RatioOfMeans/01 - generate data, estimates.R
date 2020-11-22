
require(dplyr)

dir.create("~/Temp/ratio/data", recursive = T)

set.seed(1234)

qs <- c(5,10,20,50,100)
ms <- c(5,10,20,50,100)
qs <- sort(unique(c(qs, qs^2)))
trials <- 10000
alphas <- c(0.1,1,10)

configs <- expand.grid(
  q=qs,
  m=ms,
  a=alphas
) %>%
  mutate(
    n=q*m,
    memsize=q*trials
  )

configs <- lapply(1:nrow(configs), function(i) {
  need <- ceiling(configs$memsize[i]/min(configs$memsize))
  cbind(configs[rep(i, need),], index=1:need)
}) %>%
  data.table::rbindlist() %>%
  mutate(trials=ceiling(trials*min(memsize)/memsize)) %>%
  as.data.frame()

# configs %>%
#   group_by(q, m, a) %>%
#   summarise(x=sum(trials), m=min(trials), M=max(trials)) %>%
#   as.data.frame() #%>%
  # select(x) %>%
  # table()

# end

# for balancing workload
configs <- configs[sample(1:nrow(configs), nrow(configs)),]

pbapply::pbapply(configs, 1, function(v) {
  q <- as.numeric(v["q"])
  m <- as.numeric(v["m"])
  a <- as.numeric(v["a"])
  index <- as.numeric(v["index"])
  trials <- as.numeric(v["trials"])
  tau <- 1
  
  S1 <- S2 <- S11 <- S12 <- S22 <- matrix(ncol = trials, nrow = q)
  for (j in 1:q) {
    X1 <- matrix(rgamma(m*trials, a), nrow = m, ncol = trials)
    X2 <- matrix(rgamma(m*trials, a), nrow = m, ncol = trials)
    S1[j,] <- colSums(X1)
    S2[j,] <- colSums(X2)
    S11[j,] <- colSums(X1^2)
    S12[j,] <- colSums(X1*X2)
    S22[j,] <- colSums(X2^2)
  }
  S2 <- tau*S2
  S12 <- tau*S12
  S22 <- tau^2*S22
  
  pe <- colSums(S2/S1)/q
  nc <- (pe + colSums(S12/S1^2)/q)/(1 + colSums(S11/S1^2)/q)
  ncind <- pe/(1 + colSums(S11/S1^2 - 1/m)/q)

  v <- cbind(
    "0"=colSums(S22/S1^2),
    "1"=-2*colSums(S12/S1^2),
    "2"=colSums(S11/S1^2)
  )/q^2

  res <- cbind(
    pe=pe,
    nc=nc,
    ncind=ncind,
    spe=sqrt(v[,"0"]*pe^2 + v[,"1"]*pe + v[,"2"]),
    snc=sqrt(v[,"0"]*nc^2 + v[,"1"]*nc + v[,"2"]),
    sncind=sqrt(v[,"0"]*ncind^2 + v[,"1"]*ncind + v[,"2"])
  )

  save(q, m, a, trials, tau, res, file=paste0("~/Temp/ratio/data/ratioq", q, "m", m, "a", a, "take", index, ".RData"))
})
