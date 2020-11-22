
#' Example: estimation of binomial regression

dir <- "binomial"

library(pbapply)
library(brglm2nc)
library(parallel)

meths <- c("ML", "AS_mean", "AS_median", "AS_nc")
n.sims <- 1000

confs <- expand.grid(k=c(0.1,0.3,0.5), n=c(20, 50, 100))
confs$p <- confs$k*confs$n

set.seed(12345)

# beta0 = 0, betaj = 1, rescaling covariates so that V(xb)=1 E(xb)=0
Xmat <- list()
for (i in 1:nrow(confs)) {
  conf <- confs[i,]
  Xmat[[i]] <- with(conf, matrix(rnorm(p*n)/sqrt(p), ncol = p, nrow = n))
}

cl <- makeCluster(detectCores())

rbinom2 <- Vectorize(rbinom, "prob")

res <- pbsapply(X=1:(n.sims*nrow(confs)), FUN=function(i, Xmat, rbinom2, meths, confs) {
  require(brglm2mmle)
  conf <- i %% nrow(confs) + 1
  y <- rbinom2(n=1, size=1, prob=1/(1+exp(-rowSums(Xmat[[conf]]))))
  c(conf, as.vector(sapply(meths, function(m, y, covs) {
    tryCatch({
      part <- glm(y ~ ., data = covs, family = binomial(), method = "brglmFit", type=m)
      if (part$converged) {
        phi <- part$dispersion
        part <- summary(part)
        part <- part$coefficients[1:2,1:2]
        colnames(part) <- c("est", "se")
        rownames(part) <- NULL
        part <- as.data.frame(part)
        part$bias <- part$est - c(0,1)
        part$z <- part$bias/part$se
        c(phi, unlist(part))
      } else {
        rep(NA, 9)
      }
    }, error=function(e) rep(NA, 9))
  }, y=y, covs=as.data.frame(Xmat[[conf]]))))
}, Xmat=Xmat, rbinom2=rbinom2, meths=meths, confs=confs, cl=cl)

stopCluster(cl)

res <- t(res)
colnames(res) <- c("config", apply(expand.grid(rows=c("phi", "e1", "e2", "s1", "s2", "b1", "b2", "z1", "z2"), cols=sub("^AS_","",meths)), 1, function(v) paste(rev(v), collapse = ".")))
res <- as.data.frame(res)

save(meths, res, confs, file=paste0(dir, "/joint.RData"))

meths2 <- gsub("\\.", "", unique(na.omit(stringr::str_extract(colnames(res), ".*\\."))))
newcols <- c("method", unique(sub(".*\\.", "", colnames(res))))
res2 <- as.data.frame(matrix(NA, nrow = 0, ncol = length(newcols)))
colnames(res2) <- newcols
for (m in meths2) {
  aux <- res[,!grepl("\\.", colnames(res)) | grepl(paste0("^", m, "\\."), colnames(res))]
  colnames(aux) <- gsub(".*\\.", "", colnames(aux))
  aux$method <- m
  res2 <- rbind(res2, aux)
}

save(meths2, res2, file=paste0(dir, "/pooled.RData"))
