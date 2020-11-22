source("lib.R")

verbose <- F

pbapply::pblapply(
  X=list.files(path = "~/Temp/tsn/theta", full.names = T),
  FUN=function(f, verbose, models, estResModel) {
    tryCatch({
      load(f)
      load(sub("/Temp/tsn/theta/", "/Temp/tsn/data/", f))
      mod <- models[[model]]
      for (label in colnames(res)[grepl("\\.est$", colnames(res))]) {
        est <- res[[label]]
        se <- ifelse(
          is.finite(est),
          models$general$se(
            theta = est,
            z = if (grepl("C\\.", label)) Z else Z/sqrt(1 - d/ncol(Z)),
            model = mod
          ),
          NA
        )
        res[[sub("\\..*?$", ".se", label)]] <- se
      }
      save(model, n, d, theta, res, file=f)
    }, error=function(e) save(e, file = paste0(f, ".error")))
  },
  verbose=verbose,
  models=models,
  estResModel=estResModel,
  cl = parallel::detectCores()
)
