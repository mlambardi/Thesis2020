require(dplyr)
require(ggplot2)

dat <- pbapply::pblapply(list.files("~/Temp/ratio/data", full.names = T), function(f) {
  load(f)
  as.data.frame(cbind(alpha=a, tau=tau, m=m, q=q, res))
}) %>%
  data.table::rbindlist() %>%
  as.data.frame() %>%
  mutate(
    infinite=m*alpha <= 2,
    bwpe=1/(1 - 1/(m*alpha))/sqrt(2/q*m*alpha),
    bwnc=((1/(1 - 1/(m*alpha)) + 1/m)/(1 + (1 + 1/alpha)/m) - 1)/sqrt(2/q*m*alpha),
    alpha=as.factor(alpha),
    zpe=(pe - tau)/spe,
    znc=(nc - tau)/snc,
    zncind=(ncind - tau)/sncind
  )

dat %>%
  group_by(alpha, tau, m, q) %>%
  summarise(n=n()) %>%
  as.data.frame() %>%
  select(n) %>%
  table()

for (k in 1:2) {
  dat %>%
    filter(q %in% m^k) %>%
    group_by(alpha, tau, q, m)  %>%
    summarise(
      PE=mean(zpe, na.rm=T),
      NC1=mean(znc, na.rm=T),
      NC2=mean(zncind, na.rm=T),
      `$\\pm \\delta$`=qnorm(0.975)*max(sd(zpe, na.rm=T), sd(znc, na.rm=T), sd(zncind, na.rm=T))/sqrt(n())
    ) %>%
    rename(`$\\alpha$`=alpha, `$q$`=`q`, `$m$`=m) %>%
    as.data.frame() %>%
    select(-tau) %>%
    knitr::kable(
      booktabs = T, format = "latex", escape = F, digits = 3,
      caption = "Mean bias in Wald's Z for $\\interest$. Only the largest 95\\% confidence radius is reported.",
      label = paste0("esratiosbias", k)
    ) %>%
    writeLines(paste0("output/bias", k, ".tex"))
}

thr <- 3

source("../colorblindpalette.R")

for (k in 1:2) {
  vals <- sort(unique(dat$m))^k
  for (v in 1:length(levels(dat$alpha))) {
    png(paste0("output/qqplot", k, v,".png"), res=600, units = "cm", width = 20, height = 20)
    v <- levels(dat$alpha)[v]
    dat2 <- dat %>%
      filter((q %in% vals) & (alpha == v)) %>%
      mutate(q=factor(q, levels = as.character(vals)))
    print(
      dat2 %>%
        ggplot() +
        facet_grid(q ~ m, labeller = function(...) label_both(..., sep=" = ")) +
        geom_abline(intercept = 0, slope = 1, lty=2) +
        geom_qq(aes(sample=zpe, group="PE", color="PE", lty=infinite), geom="line") +
        geom_qq(aes(sample=znc, group="NC1", color="NC1", lty=infinite), geom="line") +
        geom_qq(aes(sample=zncind, group="NC2", color="NC2", lty=infinite), geom="line") +
        xlab("theoretical") +
        ylab("sample") +
        xlim(-thr, +thr) +
        ylim(-thr, +thr) +
        geom_abline(
          data = dat2 %>%
            filter(!infinite) %>%
            group_by(q, m) %>%
            summarise(bw=mean(bwpe)),
          aes(intercept=bw, slope = 1, color="PE"), lty=2
        ) +
        theme_bw() +
        coord_fixed() +
        scale_color_manual(values = cbPalette[c(3:4,2)]) +
        ggtitle(paste0("α = ", v)) +
        guides(lty=FALSE)
    )
    dev.off()
  }
}
