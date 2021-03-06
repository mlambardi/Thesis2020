
dir <- "output"

set.seed(1234)
require(brglm2nc)

p <- 2^6
g <- sqrt(5)
k <- 0.2
n <- round(p/k)
beta <- c(0,rep(1,p))

write(p, file = paste0(dir, "/p.tex"))
write(k, file = paste0(dir, "/k.tex"))
write(n, file = paste0(dir, "/n.tex"))
write(g^2, file = paste0(dir, "/g2.tex"))
write(1, file = paste0(dir, "/minblocks.tex"))
write(p+1, file = paste0(dir, "/maxblocks.tex"))

X <- cbind(1, matrix(rnorm(n*p), ncol = p))*k/sqrt(sum(beta^2))
colnames(X) <- paste0("X", 0:p)

coefs <- matrix(NA, ncol = 5, nrow = p+1)
rownames(coefs) <- colnames(X)
colnames(coefs) <- c("PE", "NC\n~MT", "NC\nfinest", "MU", "BR")

y <- Vectorize(rbinom, "prob")(n=1, size=1, prob=1/(1 + exp(-drop(X%*%beta))))
coefs[,1] <- coef(glm(y ~ X-1, family = binomial()))
for (i in 1:(p+1)) {
  writeLines(paste(i, "of", (p+1)))
  coefs[i,2] <- coef(glm(y ~ X-1, family = binomial(), method = "brglmFit", type="AS_nc", blocks=list(alone=i, rest=setdiff(1:(p+1), i))))[i]
}
coefs[,3] <- coef(glm(y ~ X-1, family = binomial(), method = "brglmFit", type="AS_nc", blocks=1:(p+1)))
coefs[,4] <- coef(glm(y ~ X-1, family = binomial(), method = "brglmFit", type="AS_median"))
coefs[,5] <- coef(glm(y ~ X-1, family = binomial(), method = "brglmFit", type="AS_mean"))
coefs <- t(coefs)
coefs <- as.data.frame(coefs)

library(reshape2)
require(ggplot2)

#ggplot needs a dataframe
data <- coefs
#id variable for position in matrix 
data$id <- 1:nrow(data) 
#reshape to long format
data <- melt(data,id.var="id")
colnames(data)[2] <- "predictor"

source("../colorblindpalette.R")

#plot
png(paste0("output/eachaloneinturn.png"), res=600, units = "cm", width = 20, height = 20)
ggplot(data, aes(x=id,y=value,group=predictor,colour=predictor)) +
  # geom_point()+
  geom_line()+
  ylab("estimated coefficient")+
  scale_x_discrete(name ="", limits=rownames(coefs)) +
  theme_bw() +
  theme(legend.position = "none") +
  scale_colour_manual(values = c("black", rep(cbPalette,1000)))
dev.off()

res <- list()
for (i in 0:log(p,2)) {
  blocks <- split(2:(p+1), rep(1:(2^i), each=p/2^i))
  blocks$intercept <- 1
  res[[length(res)+1]] <- glm(y ~ X-1, family = binomial(), method = "brglmFit", type="AS_nc", blocks=blocks)
}
res <- t(sapply(res, coef))
res <- as.data.frame(res)

data <- res
#id variable for position in matrix 
data$id <- 1:nrow(data) 
#reshape to long format
data <- melt(data,id.var="id")
colnames(data)[2] <- "predictor"

#plot
png(paste0("output/bisect.png"), res=600, units = "cm", width = 20, height = 20)
ggplot(data, aes(x=id,y=value,group=predictor,colour=predictor)) +
  # geom_point()+
  geom_line()+
  ylab("estimated coefficients")+
  scale_x_discrete(name ="number of slopes blocks", limits=as.character(2^c(0:log(p,2))))+
  scale_colour_manual(values = c("black", rep(cbPalette,1000))) +
  theme_bw() +
  theme(legend.position = "none")
dev.off()

