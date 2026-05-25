rm(list = ls())

# libraries ---------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, tidyverse, survival, ggplot2, splines,
               mstate, rjags, R2jags, future.apply)

n_cores <- detectCores() - 1  

# load  dataset 
ns = 1
params = 1
nom <- paste0("sim_n", ns, "p", params)
load(paste0(nom, ".Rdata"))  


# jags function
bnp_simmodel <- "
  model {
    for (g in 1:G) {
      B[1, g] <- b[1, g]
    }

    for (h in 2:H) {
      for (g in 1:(G-1)) {
        B[h, g] <- b[h, g] * ((1 - b[h-1, g]) * B[h-1, g] / b[h-1, g])
      }
      B[h, G] <- 1 - sum(B[h, 1:(G-1)])
    }

    for (h in 1:H) {
      for (g in 1:G) {
        b[h, g] ~ dbeta(1, M)
        BetaProb[h, g] <- B[h, g] / sum(B[, g])
      }
    }

    for (h in 1:H) {
      for (g in 1:G) {
        betah[1:p, g, h] ~ dmnorm(mub, Sigmab)
      }
    }

    for (l in 1:n) {
      for (j in 1:ID[l]) {
        m[l, j] ~ dcat(BetaProb[, transition[l, j]])
      }
    }

    for (l in 1:n) {
      for (j in 1:ID[l]) {
        is.censored[l, j] ~ dinterval(time[l, j], cen[l, j])
        time[l, j] ~ dlnorm(mu[l, j], tau2[transition[l, j]])
        mu[l, j] <- inprod(betah[, transition[l, j], m[l, j]], X[l, ])
      }
    }

    M ~ dgamma(a0, b0)
    mub ~ dmnorm(m0, S0)
    Sigmab ~ dwish(Psiinv, psi)

    for (g in 1:G) {
      sigma2[g] <- 1 / tau2[g]
      tau2[g] ~ dgamma(gamag[g], deltag[g])
    }
  }
"

# Define function for parallel execution with logging
process_dataset <- function(data, model_str, nom, log_file, 
                            n.iter = 200000, n.burnin = 100000, n.thin = 10) {
  
  tryCatch({
    start_time <- Sys.time()
    write_log(paste0("[", format(start_time, "%Y-%m-%d %H:%M:%S"), 
                     "] Starting dataset h=", data$h))
    
    # hyperparameters
    g0 <- data$N
    Sigma0 <- g0 * solve(t(data$X) %*% data$X)
    
    # data for jags model
    d.bnp <- list(
      n = data$N, time = data$time.test, cen = data$cen.test, 
      ID = data$Ni, X = data$X, is.censored = data$is.censored.test,
      transition = data$trans.test, p = ncol(data$X), H = 25,
      gamag = rep(10,3), deltag = rep(2,3), a0 = 12, b0 = 2, G = 3,
      m0 = rep(0,3), S0 = Sigma0, Psiinv = Sigma0, psi = ncol(data$X) + 2
    )
    
    p.bnp <- c("betah", "M", "sigma2", "BetaProb")
    
    # jags run
    out.sim <- jags(
      data = d.bnp, parameters = p.bnp, model = textConnection(model_str),
      n.chains = 2, n.iter = n.iter, n.thin = n.thin, n.burnin = n.burnin
    )
    
    #results
    betah2 <- out.sim$BUGSoutput$sims.list$betah
    wh2 <- out.sim$BUGSoutput$sims.list$BetaProb
    sig2 <- out.sim$BUGSoutput$sims.list$sigma2
    M2 <- out.sim$BUGSoutput$sims.list$M
    
    save(betah2, wh2, sig2, M2, file = paste0(nom, "_h", data$h, ".Rdata"))
    
    end_time <- Sys.time()
    write_log(paste0("[", format(end_time, "%Y-%m-%d %H:%M:%S"), 
                     "] Completed dataset h=", data$h, 
                     " (Elapsed Time: ", round(difftime(end_time, start_time, units = "mins"), 2), " min)"))
    
  }, error = function(e) {
    write_log(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
                     "] ERROR in dataset h=", data$h, ": ", e$message))
  })
}


plan(multisession, workers = n_cores)  

results <- future_lapply(datas, process_dataset, 
                         model_str = bnp_simmodel, 
                         nom = nom, 
                         log_file = log_file)

# reset parallel
plan(sequential)

write_log("All datasets processed.")
