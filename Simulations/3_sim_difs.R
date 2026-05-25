rm(list = ls())

# parameters and sample size evaluation -------------------------------
ns = 3
params = 4
nom <- paste0("sim_n", ns, "p", params)

# libreries ---------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, tidyverse, survival, ggplot2, splines,
               mstate, rjags, R2jags, future.apply, parallel)

# cores -------------------------------------------------------------------
n_cores <- detectCores() - 1  

# load datas --------------------------------------------------------------
load(paste0(nom, ".Rdata"))


# objects definition ------------------------------------------------------
tt.1 <- seq(0.001, 1.35, length.out = 20)

betas <- matrix(c(2,1,-1,-1,1,2,
                  2,2,2,2,2,2,
                  -0.03,-0.49,-0.01,0.15,0.002,0.04,
                  rep(0.001,6)),
                ncol = 6,nrow=4,byrow = T)

rates <- matrix(c(1,1.5,1,2,1,2.5,
                  1.5,1.5,1.5,1.5,1.5,1.5,
                  2.59 ,0.39,1.52,0.66,1.52,0.66,
                  rep(1,6)),
                ncol = 6,nrow=4,byrow = T)

# log file ----------------------------------------------------------------
log_file <- paste0(nom, "difs_log.txt")


# parallel execution -----------------------------------------
for (index in 1:100) {
  
  start_time <- Sys.time()
  write_log(paste0("[", format(start_time, "%Y-%m-%d %H:%M:%S"), 
                   "] Starting dataset h=", index))
  
  # Check if dataset file exists
  dataset_file <- paste0(nom, "_h", index, ".Rdata")
  
  load(dataset_file)
  
  
  S_pred <- function(y, xx, tr, betah2, sig2, wh2, Nmc, Hdp) {
    meanlog <- vapply(1:Hdp, function(jdp) {
      betah2[,, tr, jdp] %*% xx  
    }, numeric(Nmc))
    
    sdlog <- sqrt(sig2[, tr])  
    
    plnorm_terms <- vapply(1:Hdp, function(jdp) {
      plnorm(y, meanlog[, jdp], sdlog)
    }, numeric(Nmc))
    
    Nmc <- dim(betah2)[1]
    Hdp <- dim(betah2)[4]
    tra <- dim(betah2)[3]
    p <- dim(betah2)[2]
    
    weighted_sum <- rowSums(wh2[, , tr] * plnorm_terms)  
    
    return(1 - weighted_sum)  
  }
  
  
  S_pred_mean <- function(y, xx, tr, betah2, sig2, wh2, Nmc, Hdp, cores = detectCores() - 1) {
    
    eval.pred <- mclapply(1:length(y), function(h) {
      S_pred(y[h], xx, tr, betah2, sig2, wh2, Nmc, Hdp)
    }, mc.cores = cores)
    
    eval.pred.mean <- sapply(eval.pred, mean)  
    return(eval.pred.mean)
  }
  
  Nmc <- dim(betah2)[1]
  Hdp <- dim(betah2)[4]
  tra <- dim(betah2)[3]
  p <- dim(betah2)[2]
  
  mssimdat <- datas[[index]]$mssimdat
  mssimdat$trans <- factor(mssimdat$trans)
  
  
  # cox model
  c2 <- coxph(Surv(time, status) ~ strata(trans) + x1.1 + x1.2 + x1.3 + x2.1 + x2.2 + x2.3, data = mssimdat)
  
  # hazards for all transitions in the dataset
  Hazc2.c <- basehaz(c2, centered = FALSE)
  all_transitions <- unique(mssimdat$trans)
  Hazc2_strata <- lapply(all_transitions, function(tr) {
    Hazc2.c[Hazc2.c$strata == tr, ]
  })
  names(Hazc2_strata) <- all_transitions
  
  library(parallel)
  
  # process for one subject
  process_subject <- function(h) {
    nd <- mssimdat[mssimdat$id == h, ]
    Xnew <- c(1, unique(nd$x1), unique(nd$x2))
    trans.nums <- as.numeric(nd$trans)
    
    difs_subject <- matrix(0, nrow = length(trans.nums), ncol = 4)
    
    for (i in seq_along(trans.nums)) {
      tr <- trans.nums[i]
      
      haz <- Hazc2_strata[[as.character(tr)]]
      baseline_hazard <- haz$hazard
      
      lin_pred <- sum(coef(c2) * nd[i, names(coef(c2))])
      survival_prob <- exp(-baseline_hazard * exp(lin_pred))
      
      tt.1 <- haz$time
      
      bnp_prob <- S_pred_mean(tt.1, xx = Xnew, tr = tr, betah2, sig2, wh2, Nmc, Hdp)
      
      real <- exp(-(tt.1^rates[params, 2 * tr]) * rates[params, (2 * tr - 1)] * 
                    exp(betas[params, (2 * tr - 1)] * Xnew[2] + betas[params, 2 * tr] * Xnew[3]))
      
      difs_subject[i, ] <- c(
        sum(abs(real - survival_prob)),
        sum(abs(real - bnp_prob)),
        sum((real - survival_prob)^2),
        sum((real - bnp_prob)^2)
      )
    }
    #sum transitions for each subject
    return(colSums(difs_subject)) 
  }
  
  # parallel for all the datasets
  cl <- makeCluster(detectCores() - 1)
  clusterExport(cl, c("process_subject", "S_pred_mean", "betah2", "sig2", "wh2", "Nmc", "Hdp", "S_pred",
                      "Hazc2_strata", "c2", "rates", "params", "betas", "mssimdat","mclapply","detectCores"))
  
  difs.s <- parLapply(cl, 1:length(unique(mssimdat$id)), process_subject)
  stopCluster(cl)
  
  res <- sqrt(colSums(do.call(rbind, difs.s), na.rm = TRUE))
  save(res, file = paste0(nom, "_difs_h", index,"_new.Rdata"))
}
