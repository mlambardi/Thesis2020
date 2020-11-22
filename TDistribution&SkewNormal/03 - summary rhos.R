
require(dplyr)
require(ggplot2)

dat <- pbapply::pblapply(list.files("~/Temp/tsn/rho", full.names = T), function(f) {
  load(f)
  rho <- as.data.frame(rho)
  cbind(rho, d=d, n=n, model=model, theta=theta)
}) %>%
  data.table::rbindlist() %>%
  arrange(model, theta) %>%
  mutate(label=factor(paste0(model, ifelse(theta < 10, "0", ""), theta))) %>%
  as.data.frame()

ks <- sort(as.numeric(unique(sub("^rho(\\d+).*?$", "\\1", colnames(dat)[grepl("^rho\\d+", colnames(dat))]))))
mets <- c("pe", "naive", "nc2", "nc3onestep", "nc3")

for (m in mets) {
  r3 <- dat[[paste0("rho3", m)]]
  r4 <- dat[[paste0("rho4", m)]]
  for (k in 3:4) {
    rk <- dat[[paste0("rho", k, m)]]
    rkm1 <- if (k==3) 1 else dat[[paste0("rho", k-1, m)]]
    rkp1 <- dat[[paste0("rho", k+1, m)]]
    rkp2 <- dat[[paste0("rho", k+2, m)]]
    r2k <- dat[[paste0("rho", 2*k, m)]]
    dat[[paste0("serho", k, m)]] <- sqrt((k^2*rkm1^2 + r3*k^2*rkm1*rk - 2*k*rkm1*rkp1 + k^2*rk^2*(r4-1)/4 - k*rk*(rkp2 - rk) + r2k - rk^2)/dat$n)
  }
}

summary(dat[,grepl("^se", colnames(dat))])
summary(dat[,grepl("^rho4", colnames(dat))] > 1)

for (k in 3:4) {
  for (m1 in 2:length(mets)) {
    for (m2 in mets[(m1-1):1]) {
      dat[[paste0("serho", k, mets[m1])]] <- ifelse(
        is.na(dat[[paste0("serho", k, mets[m1])]]),
        dat[[paste0("serho", k, m2)]],
        dat[[paste0("serho", k, mets[m1])]]
      )
    }
  }
}

source("lib.R")

for (k in 3:4) {
  dat[[paste0("rho", k, "true")]] <- ifelse(
    dat$model=="t",
    models$t$rho(k = k, theta = dat$theta),
    models$sn$rho(k = k, theta = dat$theta)
  )
  for (m in mets) {
    rk <- dat[[paste0("rho", k, m)]]
    theta <- dat[[paste0("rho", k, "true")]]
    se <- dat[[paste0("serho", k, m)]]
    dat[[paste0("wrho", k, m)]] <- (rk - theta)/se
    dat[[paste0("drho", k, m)]] <- (rk - theta)*sqrt(dat$n)
  }
}

source("../colorblindpalette.R")

thr <- 3

for (lab in levels(dat$label)) {
  for (k in 3:4) {
    writeLines(paste0(lab, " ", k))
    model <- switch(dat$model[dat$label==lab][1], "t"="t-distribution", "skew-normal distribution")
    theta <- dat$theta[dat$label==lab][1]
    rho <- dat[dat$label==lab,paste0("rho",k,"true")][1]
    png(paste0("output/rho", k, which(levels(dat$label)==lab),".png"), res=600, units = "cm", width = 20, height = 20)
    print(
      dat %>%
        filter(label==lab) %>%
        ggplot() +
        facet_grid(d ~ n/d, labeller = function(...) label_both(..., sep=" = ")) +
        geom_abline(intercept = 0, slope = 1, lty=2) +
        xlab("theoretical") +
        ylab("sample") +
        xlim(-thr, +thr) +
        ylim(-thr, +thr) +
        theme_bw() +
        theme(legend.position="top") +
        coord_fixed() +
        scale_color_manual(values = cbPalette) +
        ggtitle(paste0("ρ", k, "=", signif(rho,3)), subtitle=paste0(model, ", τ=", theta)) +
        geom_qq(aes_string(sample=paste0("wrho",k,"pe"), group="'PE'", color="'PE'"), geom="line") +
        geom_qq(aes_string(sample=paste0("wrho",k,"naive"), group="'naive'", color="'naive'"), geom="line") +
        geom_qq(aes_string(sample=paste0("wrho",k,"nc2"), group="'NC2'", color="'NC2'"), geom="line") +
        geom_qq(aes_string(sample=paste0("wrho",k,"nc3"), group="'NC3'", color="'NC3'"), geom="line")
    )
    dev.off()
  }
}

stdbias <- dat %>%
  group_by(model, theta) %>%
  summarise(
    rho3pe=mean(rho3pe-rho3true, na.rm = T)/sd(rho3pe-rho3true, na.rm = T),
    rho3nc2=mean(rho3nc2-rho3true, na.rm = T)/sd(rho3nc2-rho3true, na.rm = T),
    rho3nc3=mean(rho3nc3-rho3true, na.rm = T)/sd(rho3nc3-rho3true, na.rm = T),
    rho4pe=mean(rho4pe-rho4true, na.rm = T)/sd(rho4pe-rho4true, na.rm = T),
    rho4nc2=mean(rho4nc2-rho4true, na.rm = T)/sd(rho4nc2-rho4true, na.rm = T),
    rho4nc3=mean(rho4nc3-rho4true, na.rm = T)/sd(rho4nc3-rho4true, na.rm = T)
  )

stdbias[,-1] <- apply(stdbias[,-1], 2, function(x)as.character(signif(x,digits=2)))
