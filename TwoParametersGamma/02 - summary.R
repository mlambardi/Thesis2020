require(tidyr)
require(dplyr)
require(ggplot2)
require(ggrepel)
require(knitr)
options(knitr.table.format = "latex")

res <- pbapply::pblapply(list.files(path = "~/Temp/gamma/data", full.names = T), function(f) {
  load(f)
  cbind(n=n, a=a, l=l, trials=trials, res)
}) %>%
  data.table::rbindlist() %>%
  as.data.frame()

# res2 <- res %>%
#   group_by(α, n) %>%
#   summarise(
#     xα=c(
#       qnorm(mean(wαpe < 3/5, na.rm=T)),
#       qnorm(mean(wαnc < 1/5, na.rm=T)),
#       qnorm(mean(wαbr < -1/5, na.rm=T)),
#       qnorm(mean(wαbrlog < -3/5, na.rm=T))
#     ),
#     xλ=c(
#       qnorm(mean(wλpe < 3/5, na.rm=T)),
#       qnorm(mean(wλnc < 1/5, na.rm=T)),
#       qnorm(mean(wλbr < -1/5, na.rm=T)),
#       qnorm(mean(wλbrlog < -3/5, na.rm=T))
#     ),
#     y=c(
#       3/5, 1/5, -1/5, -3/5
#     ),
#     label=c(
#       "PE",
#       "NC",
#       "BR",
#       "BR-log"
#     )
#   ) %>%
#   as.data.frame()

meths <- c("pe", "nc", "br", "brlog")
for (m in meths) {
  criterion <- paste0("scorecrit", m)
  alpha.est <- paste0("a", m)
  lambda.est <- paste0("l", m)
  alpha.se <- paste0("sea", m)
  lambda.se <- paste0("sel", m)
  alpha.wald <- paste0("wa", m)
  lambda.wald <- paste0("wl", m)
  res[[alpha.est]] <- ifelse(is.finite(res[[alpha.est]]), res[[alpha.est]], NA)
  res[[lambda.est]] <- ifelse(is.finite(res[[lambda.est]]), res[[lambda.est]], NA)
  conv <- res[[criterion]]
  conv <- ifelse(is.finite(conv), abs(conv) < 1e-7, F)
  res[[alpha.est]][!conv] <- NA
  res[[lambda.est]][!conv] <- NA
  res[[alpha.wald]] <- res[[alpha.est]]*log(res[[alpha.est]]/res$a)/res[[alpha.se]]
  res[[lambda.wald]] <- res[[lambda.est]]*log(res[[lambda.est]]/res$l)/res[[lambda.se]]
}
res <- res[,c("a", "l", "n", apply(expand_grid(c("a", "l", "sea", "sel", "wa", "wl"), meths), 1, paste, collapse=""))]
colnames(res) <- sub("λog", "log", sub("a", "α", sub("l", "λ", colnames(res))))

rm(list = setdiff(ls(), "res"))
gc()

source("../colorblindpalette.R")

res %>%
  group_by(α, n) %>%
  summarize(
    αpe=100*mean(!is.finite(αpe)),
    αnc=100*mean(!is.finite(αnc)),
    αbr=100*mean(!is.finite(αbr)),
    αbrlog=100*mean(!is.finite(αbrlog)),
  	trials=n()
  ) %>%
  filter(pmax(αpe, αnc, αbr, αbrlog) > 0) %>%
  rename(`$\\alpha$`=α, PE=αpe, NC=αnc, BR=αbr, `BR-log`=αbrlog) %>%
  select(-PE, -NC, -`BR-log`, -trials) %>%
  as.data.frame() %>%
  knitr::kable(
    booktabs = T, format = "latex", escape = F, digits = 3,
    caption = "Non-convergence rate ($\\%$) in 100000 simulations, gamma case.",
    label = "esgammaconv"
  ) %>%
  writeLines("output/convergence.tex")

out <- res %>%
  group_by(α, n) %>%
  summarise(
    `$\\alpha$ : PE`=mean(wαpe, na.rm=T),
    NC=mean(wαnc, na.rm=T),
    BR=mean(wαbr, na.rm=T),
    `BR-log`=mean(wαbrlog, na.rm=T),
    aux1=max(sd(wαpe, na.rm=T), sd(wαnc, na.rm=T), sd(wαbr, na.rm=T), sd(wαbrlog, na.rm=T))/sqrt(n()),
    `$\\lambda$ : PE`:=mean(wλpe, na.rm=T),
    ` NC`=mean(wλnc, na.rm=T),
    ` BR`=mean(wλbr, na.rm=T),
    ` BR-log`=mean(wλbrlog, na.rm=T),
    aux2=max(sd(wλpe, na.rm=T), sd(wλnc, na.rm=T), sd(wλbr, na.rm=T), sd(wλbrlog, na.rm=T))/sqrt(n())
    # , trials=n()
  ) %>%
  rename(`$\\alpha$`=α)

out %>%
  select(-aux1, -aux2) %>%
  as.data.frame() %>%
  knitr::kable(
    booktabs = T, format = "latex", escape = F, digits = 3,
    caption = paste0("Estimated mean bias of Wald's Z statistic for $\\log\\alpha$ and $\\log\\lambda$, gamma case. All standard errors are all less or equal than $", signif(max(c(out$aux1, out$aux2)), 1), "$."),
    label = "esgammabias"
  ) %>%
  as.character() %>%
  gsub(pattern = "\\\\addlinespace", replacement = "") %>%
  sub(pattern = "rrrrrrrrrr", replacement = "rr|rrrr|rrrr") %>%
  writeLines("output/bias.tex")

gc()

thr <- 3

# jpeg("output/loglambda-qqplotbw.jpg", res=300, units = "cm", width = 20, height = 20)
# res %>%
#   ggplot() +
#   facet_grid(α ~ n, labeller = function(...) label_both(..., sep=" = ")) +
#   geom_abline(intercept = 0, slope = 1, lty=2) +
#   geom_qq(aes(sample=wλpe, group="PE"), geom="line") +
#   geom_qq(aes(sample=wλnc, group="NC"), geom="line") +
#   geom_qq(aes(sample=wλbr, group="BR"), geom="line") +
#   geom_qq(aes(sample=wλbrlog, group="BR-log"), geom="line") +
#   xlim(-thr,+thr) +
#   ylim(-thr,+thr) +
#   theme(legend.position="top") +
#   coord_fixed() +
#   theme_bw() +
#   geom_text_repel(data = res2, aes(x=xλ, y=y, label=label))
# dev.off()

slides <- F

jpeg(paste0("output/loglambda-qqplot", if (slides) "simple" else "",".jpg"), res=600, units = "cm", width = 20, height = 20-6*slides)
res %>%
  filter(!slides | (α %in% c(0.1, 1, 10))) %>%
  ggplot() +
  facet_grid(α ~ n, labeller = function(...) label_both(..., sep=" = ")) +
  geom_abline(intercept = 0, slope = 1, lty=2) +
  geom_qq(aes(sample=wλpe, group="PE", color="PE"), lty=2, geom="line") +
  geom_qq(aes(sample=wλnc, group="NC", color="NC"), lty=2, geom="line") +
  geom_qq(aes(sample=wλbr, group="BR", color="BR"), lty=2, geom="line") +
  geom_qq(aes(sample=wλbrlog, group="BR-log", color="BR-log"), lty=2, geom="line") +
  xlim(-thr,+thr) +
  ylim(-thr,+thr) +
  theme(legend.position="top") +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = cbPalette)
dev.off()

jpeg(paste0("output/logalpha-qqplot", if (slides) "simple" else "", ".jpg"), res=600, units = "cm", width = 20, height = 20-6*slides)
res %>%
  filter(!slides | (α %in% c(0.1, 1, 10))) %>%
  ggplot() +
  facet_grid(α ~ n, labeller = function(...) label_both(..., sep=" = ")) +
  geom_abline(intercept = 0, slope = 1, lty=2) +
  geom_qq(aes(sample=wαpe, group="PE", color="PE"), lty=2, geom="line") +
  geom_qq(aes(sample=wαnc, group="NC", color="NC"), lty=2, geom="line") +
  geom_qq(aes(sample=wαbr, group="BR", color="BR"), lty=2, geom="line") +
  geom_qq(aes(sample=wαbrlog, group="BR-log", color="BR-log"), lty=2, geom="line") +
  xlim(-thr,+thr) +
  ylim(-thr,+thr) +
  theme(legend.position="top") +
  coord_fixed() +
  theme_bw() +
  scale_color_manual(values = cbPalette)
dev.off()
