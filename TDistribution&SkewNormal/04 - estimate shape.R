source("lib.R")

verbose <- F

if (!dir.exists("~/Temp/tsn/theta")) {
  dir.create("~/Temp/tsn/theta", recursive = T)
}

fs.needed <- list.files(path = "~/Temp/tsn/data", full.names = T)
fs.avail <- sub("/Temp/tsn/theta/", "/Temp/tsn/data/", list.files(path = "~/Temp/tsn/theta", full.names = T))
pbapply::pblapply(
  X=setdiff(fs.needed, fs.avail),
  FUN=function(f, verbose, models, estResModel) {
    f2 <- sub("/Temp/tsn/data/", "/Temp/tsn/theta/", f)
    tryCatch({
      load(f)
      res <- as.data.frame(matrix(NA, ncol = 0, nrow = nrow(Z)))
      for (enhancesigma in c(F, T)) {
        for (corrbeta in c(F, T)) {
          for (corrsigma in (if (enhancesigma & corrbeta) c(F, T) else (enhancesigma | corrbeta))) {
            label <- paste0(ifelse(corrbeta,"B","b"), ifelse(corrsigma,"S","s"), ifelse(enhancesigma,"C","c"))
            
            if (verbose) {
              writeLines(paste0("correcting for ", paste(c(if (corrbeta) "beta", if (corrsigma) "sigma", if (!corrbeta & !corrsigma) "none (classical estimation)"), collapse = " & ")))
            }
            
            res[[paste0(label, ".est")]] <- estResModel(
              Z = Z,
              d = d,
              mod = models[[model]],
              corrbeta = corrbeta,
              corrsigma = corrsigma,
              enhancesigma = enhancesigma,
              verbose = F
            )
          }
        }
      }
      save(model, n, d, theta, res, file=f2)
    }, error=function(e) save(e, file = paste0(f2, ".error")))
  },
  verbose=verbose,
  models=models,
  estResModel=estResModel,
  cl = parallel::detectCores()
)
