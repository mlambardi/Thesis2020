  
#' Example: estimation of binomial regression - analysis of the results

rm(list=ls())
dir <- "output"

source("../colorblindpalette.R")

for (f in list.files(dir, pattern = "RData$", full.names = T)) {
  load(f)
}

library(ggplot2)
library(dplyr)

thr <- 100

aux <- grepl("\\.z\\d+$", colnames(res))
summary(res[,aux])
res[,aux][abs(res[,aux]) > thr] <- NA
summary(is.na(res[,aux]))
aux <- c("z1","z2")
res2[,aux][abs(res2[,aux]) > thr] <- NA
summary(is.na(res2[,aux]))

res$n <- confs$n[res$config]
res$k <- confs$k[res$config]
res2$n <- confs$n[res2$config]
res2$k <- confs$k[res2$config]

sz <- 0.3

require(xtable)
out <- res %>%
  group_by("p/n"=k, "n"=n) %>%
  summarise(
    PE=100*mean(is.na(ML.z1)),
    NC=100*mean(is.na(modmle.z1)),
    MU=100*mean(is.na(median.z1)),
    BR=100*mean(is.na(mean.z1)),
    sims=n()
  )
out <- xtable(out %>% select(-sims),
              digits = c(1,1,0,rep(1,4)),
              label=paste0("tab:binomialmoderatemissing"),
              caption = paste0("Non-convergence rate (\\%) in ", out$sims[1]," simulations."),
              type = "latex")
print(out, file = paste0(dir, "/missing.tex"), include.rownames=FALSE)

out <- res %>%
  group_by("p/n"=k, "n"=n) %>%
  summarise(
    PE=mean(ML.b2, na.rm=T),
    NC=mean(modmle.b2, na.rm=T),
    MU=mean(median.b2, na.rm=T),
    BR=mean(mean.b2, na.rm=T),
    sPE=sd(ML.b2, na.rm=T)/sqrt(n()),
    sNC=sd(modmle.b2, na.rm=T)/sqrt(n()),
    sMU=sd(median.b2, na.rm=T)/sqrt(n()),
    sBR=sd(mean.b2, na.rm=T)/sqrt(n()),
    sims=n()
  ) %>%
  mutate(`max. SE`=pmax(sPE, sNC, sMU, sBR)) %>%
  select(-sPE, -sNC, -sMU, -sBR)
out <- xtable(out %>% select(-sims),
              digits = c(1,1,0,rep(1,4)),
              label=paste0("tab:binomialmoderatemissing"),
              caption = paste0("Mean bias in finite estimates of first slope."),
              type = "latex")
print(out, file = paste0(dir, "/bias.tex"), include.rownames=FALSE)

# res$n <- factor(res$n)
# res$k <- factor(res$k)

png(paste0(dir, "/waldb0.png"), res=300, width=15, height =15, units = "cm")
ggplot(res) +
  facet_grid(k ~ n, labeller = function(...) label_both(..., sep=" = ")) +
  xlab("theoretical") +
  ylab("sample") +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  stat_qq(aes(sample=ML.z1  , color="PE"), geom = "line") +
  stat_qq(aes(sample=mean.z1  , color="BR"), geom = "line") +
  stat_qq(aes(sample=median.z2  , color="MU"), geom = "line") +
  stat_qq(aes(sample=modmle.z1, color="NC" ), geom = "line") +
  theme_bw() +
  theme(legend.position = "top") +
  scale_color_manual(values = cbPalette)
dev.off()

png(paste0(dir, "/waldb1.png"), res=300, width=15, height =15, units = "cm")
ggplot(res) +
  facet_grid(k ~ n, labeller = function(...) label_both(..., sep=" = ")) +
  xlab("theoretical") +
  ylab("sample") +
  geom_abline(intercept = 0, slope = 1) +
  coord_fixed() +
  stat_qq(aes(sample=ML.z2  , color="PE"), geom = "line") +
  stat_qq(aes(sample=mean.z2  , color="BR"), geom = "line") +
  stat_qq(aes(sample=median.z2  , color="MU"), geom = "line") +
  stat_qq(aes(sample=modmle.z2, color="NC" ), geom = "line") +
  theme_bw() +
  theme(legend.position = "top") +
  scale_color_manual(values = cbPalette)
dev.off()

ggplot(res, aes(x=ML.z2, xend=ML.z2)) + facet_grid(k ~ n) + xlab("n") + ylab("p/n") + ggtitle("Wald Z") +
  # geom_point(aes(y=median.z2-ML.z2, color="MU-PE" )) +
  # geom_point(aes(y=modmle.z2-ML.z2, color="NC-PE")) +
  # geom_point(aes(y=mean.z2-ML.z2, color="BR-PE")) +
  # geom_abline(intercept = 0, slope = 1) +
  geom_segment(aes(y=pmin(median.z2, modmle.z2, mean.z2)-ML.z2, yend=pmax(median.z2, modmle.z2, mean.z2)-ML.z2)) +
  coord_fixed() +
  scale_color_manual(values = cbPalette[-1])

res %>%
  group_by(k, n, meth) %>%
  summarise()
  facet_grid(k ~ n, labeller = function(...) label_both(..., sep=" = ")) +
  xlab("theoretical") +
  ylab("sample") +
  geom_abline(intercept = 0, slope = 1) +
  geom_boxplot(aes(x=sqrt(n)*b2, y=method)) +
  theme_bw() +
  theme(legend.position = "top") +
  scale_color_manual(values = cbPalette)

ggplot(res2$method, aes(x=ML.b2)) + facet_grid(k ~ n) + xlab("n") + ylab("p/n") + ggtitle("Wald Z") +
  # geom_point(aes(y=median.z2-ML.z2, color="MU-PE" )) +
  # geom_point(aes(y=modmle.z2-ML.z2, color="NC-PE")) +
  # geom_point(aes(y=mean.z2-ML.z2, color="BR-PE")) +
  # geom_abline(intercept = 0, slope = 1) +
  geom_segment(aes(y=pmin(median.z2, modmle.z2, mean.z2)-ML.z2, yend=pmax(median.z2, modmle.z2, mean.z2)-ML.z2)) +
  coord_fixed() +
  scale_color_manual(values = cbPalette[-1])

ggplot(res2, aes(group=1)) + facet_grid(k ~ n) + xlab("n") + ylab("p/n") + geom_boxplot(aes(y=b1, fill=method))
