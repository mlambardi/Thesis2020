require(ggplot2)
require(dplyr)
require(scales) # see trans_new

dat <- pbapply::pblapply(list.files("~/Temp/tsn/theta", pattern = "\\.RData$", full.names = T), function(f) {
  load(f)
  cbind(as.data.frame(res), d=d, n=n, model=model, theta=theta)
}) %>%
  data.table::rbindlist() %>%
  filter(model=="t")%>%
  rename(PE=bsc.est, MT=BSc.est, NC1=bSC.est, NC2=BsC.est, NC3=BSC.est) %>%
  mutate(config=theta) %>%
  as.data.frame()

out <- dat %>%
  group_by(theta, n/d, d) %>%
  summarise(
    "Inf : PE"=100*mean(!is.finite(PE)),
    "MT "=100*mean(!is.finite(MT)),
    "NC1 "=100*mean(!is.finite(NC1)),
    "NC2 "=100*mean(!is.finite(NC2)),
    "NC3 "=100*mean(!is.finite(NC3)),
    "Q2 : PE"=median(PE, 0.5),
    " MT"=median(MT, 0.5),
    " NC1"=median(NC1, 0.5),
    " NC2"=median(NC2, 0.5),
    " NC3"=median(NC3, 0.5),
  ) %>%
  rename("$\\tau$"=theta) %>%
  as.data.frame() %>%
  knitr::kable(
    booktabs = T, format = "latex", escape = F, digits = 1,
    caption = "Non-convergence rate (\\%) and median estimate for $\\interest$",
    label = "summarytmod"
  )

gsub("\\\\addlinespace", "", as.character(out)) %>%
  writeLines(paste0("output/tmod.tex"))

source("../colorblindpalette.R")

type <- geom_boxplot
# type <- geom_violin

for (v in unique(dat$config)) {
  theta <- as.numeric(sub("^.*?(\\d+)$", "\\1", v))
  print(v)
  png(paste0("output/box", theta, ".png"), width = 25, height = 20, units = "cm", res = 100)
  print(
    dat %>%
      filter(config==v) %>%
      mutate() %>%
      ggplot() +
      facet_grid(paste("d =", d) ~ paste("n/d =", n/d)) +
      geom_hline(yintercept = theta) +
      type(aes(y=PE, x="PE", color="PE")) +
      type(aes(y=MT, x="MT", color="MT")) +
      type(aes(y=NC3, x="NC1", color="NC1")) +
      type(aes(y=NC3, x="NC2", color="NC2")) +
      type(aes(y=NC3, x="NC3", color="NC3")) +
      scale_y_continuous(
        breaks = c(1, theta, 100),
        limits=c(2, 120)
      ) +
      scale_x_discrete(breaks=c()) +
      xlab("") +
      ylab("estimates") +
      # coord_trans(y=tn1) +
      theme_bw() +
      theme(legend.position = "top") +
      ggtitle(paste0("τ = ", v)) +
      scale_color_manual(values = cbPalette[-1])
  )
  dev.off()
}

thr <- 150
type <- geom_smooth
for (v in unique(dat$config)) {
  theta <- as.numeric(sub("^.*?(\\d+)$", "\\1", v))
  print(v)
  png(paste0("output/shrink", theta, ".png"), width = 25, height = 20, units = "cm", res = 100)
  print(
    dat %>%
      mutate(
        PE=ifelse(PE > thr, NA, PE),
        MT=ifelse(MT > thr, NA, MT),
        NC1=ifelse(NC1 > thr, NA, NC1),
        NC2=ifelse(NC2 > thr, NA, NC2),
        NC3=ifelse(NC3 > thr, NA, NC3)
      ) %>%
      filter(config==v) %>%
      ggplot() +
      facet_grid(paste("d =", d) ~ paste("n/d =", n/d)) +
      coord_fixed(xlim = c(0,thr), ylim=c(0,thr)) +
      geom_hline(yintercept = theta) +
      geom_vline(xintercept = theta) +
      geom_abline(intercept = 0, slope = 1, lty=2) +
      type(aes(x=PE, y=MT, group="MT", color="MT"), se=F) +
      type(aes(x=PE, y=NC1, group="NC1", color="NC1"), se=F) +
      type(aes(x=PE, y=NC2, group="NC2", color="NC2"), se=F) +
      type(aes(x=PE, y=NC3, group="NC3", color="NC3"), se=F) +
      xlab("classical") +
      ylab("enhanced") +
      theme_bw() +
      theme(legend.position = "top") +
      ggtitle(paste0("τ = ", v)) +
      scale_color_manual(values = cbPalette[-1])
  )
  dev.off()
}

type <- geom_point
png(paste0("output/shrink", theta, "scatter.png"), width = 25, height = 20, units = "cm", res = 100)
print(
  dat %>%
    mutate(
      PE=ifelse(PE > thr, NA, PE),
      MT=ifelse(MT > thr, NA, MT),
      NC1=ifelse(NC1 > thr, NA, NC1),
      NC2=ifelse(NC2 > thr, NA, NC2),
      NC3=ifelse(NC3 > thr, NA, NC3)
    ) %>%
    filter(config==v) %>%
    ggplot() +
    facet_grid(paste("d =", d) ~ paste("n/d =", n/d)) +
    coord_fixed(xlim = c(0,thr), ylim=c(0,thr)) +
    geom_hline(yintercept = theta) +
    geom_vline(xintercept = theta) +
    geom_abline(intercept = 0, slope = 1, lty=2) +
    type(aes(x=PE, y=MT, group="MT", color="MT"), size=0.1) +
    type(aes(x=PE, y=NC1, group="NC1", color="NC1"), size=0.1) +
    type(aes(x=PE, y=NC2, group="NC2", color="NC2"), size=0.1) +
    type(aes(x=PE, y=NC3, group="NC3", color="NC3"), size=0.1) +
    xlab("classical") +
    ylab("enhanced") +
    theme_bw() +
    theme(legend.position = "top") +
    ggtitle(paste0("τ = ", v)) +
    scale_color_manual(values = cbPalette[-1])
)
dev.off()
