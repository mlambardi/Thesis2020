
require(brglm2nc) # for br methods
require(survival) # for conditional likelihood
require(datasets) # for infertility data
require(dplyr)

print(infert)

infert$stratum <- factor(infert$stratum)
infert$spontaneous <- factor(infert$spontaneous)
infert$induced <- factor(infert$induced)

N <- length(levels(infert$stratum))

nuis <- list(
  "1"=list(1:N),
  "N"=lapply(1:N, identity)
)
interest <- c(1,2,4)

res <- list()

res[[1]] <- glm(case ~ stratum - 1 + spontaneous + induced, data = infert, family = binomial())
res[[1]]$type <- "ML"

res[[2]] <- clogit(case ~ spontaneous + induced + strata(stratum), data = infert)
res[[2]]$type <- "CL"

res[[3]] <- glm(case ~ stratum - 1 + spontaneous + induced, data = infert, family = binomial(), method = "brglmFit", type="AS_mean")
res[[3]]$type <- "Firth"

res[[4]] <- glm(case ~ stratum - 1 + spontaneous + induced, data = infert, family = binomial(), method = "brglmFit", type="AS_median")
res[[4]]$type <- "MU"

for (lambda in names(nuis)) {
  print(lambda)
  for (psi in interest) {
    print(psi)
    scheme <- nuis[[lambda]]
    if (psi==1) {
      scheme$psis <- N+1:4
    } else {
      if (psi==2) {
        scheme$psi12 <- N+1:2
        scheme$psi34 <- N+3:4
      } else {
        # psi==4
        scheme$psi1 <- N+1
        scheme$psi2 <- N+2
        scheme$psi3 <- N+3
        scheme$psi4 <- N+4
      }
    }
    i <- length(res)+1
    res[[i]] <- glm(case ~ stratum - 1 + spontaneous + induced, data = infert, family = binomial(), method = "brglmFit", type="AS_nc", blocks=scheme)
    res[[i]]$type <- "MMLE"
    res[[i]]$lambdas <- lambda
    res[[i]]$psis    <- psi
  }
}

res <- lapply(res, function(r) {
  r <- summary(r)
  r$coefficients <- r$coefficients[c("spontaneous1", "spontaneous2", "induced1", "induced2"),]
})

names(res) <- c("PE", "CL", "BR", "MU", apply(expand.grid("NC", interest, names(nuis)), 1, paste, collapse="-"))

colnames(res$CL)[c(1, 3)] <- c("Estimate", "Std. Error")

res <- t(sapply(res, function(r) apply(signif(as.matrix(r[,c("Estimate", "Std. Error")]),4), 1, function(v) paste0(v[1], " (", v[2], ")"))))

knitr::kable(res, format = "latex", label = "infertdata", caption = "Estimates of structural parameters for infertility data analysis under different approaches") %>%
  gsub(pattern = "\\\\hline", replacement="") %>%
  write(file = "output/infert.tex")
