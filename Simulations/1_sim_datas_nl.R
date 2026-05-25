#simulation additional scenarios
rm(list = ls())

# libraries ---------------------------------------------------------------
library(mstate)

# functions - load --------------------------------------------------------
source("0_sim_functions_v2.R")

# functions - aditional ---------------------------------------------------------------
rtrans <- function(n, x1, x2, pars) {
  # pars: type, lambda, k, beta1, beta2. type should be one of these options
  #"linear" exp(beta1*x1 + beta2*x2); "nonlinear"  exp(beta1*x1^2 + beta2*x2); "nph_logt"   exp(beta1*x1 + beta2*x2*log(t))
  u <- runif(n)
  if (pars$type == "linear") {
    eta <- pars$beta1 * x1 + pars$beta2 * x2
    t <- (-log(1 - u) / (pars$lambda * exp(eta)))^(1 / pars$k)
  } else if (pars$type == "nonlinear") {
    eta <- pars$beta1 * x1^2 + pars$beta2 * x2
    t <- (-log(1 - u) / (pars$lambda * exp(eta)))^(1 / pars$k)
  } else if (pars$type == "nph_logt") {
    alpha <- pars$k + pars$beta2 * x2
    t <- ((-log(1 - u)) * alpha / (pars$lambda * pars$k * exp(pars$beta1 * x1)))^(1 / alpha)
  } 
  return(t)
}

sim_data_nl <- function(n, trans_pars, cens.par, ...) {
  x1 <- runif(n)
  x2 <- rbinom(n, size = 1, prob = 0.5)
  # trans times
  timesick <- rtrans(n, x1, x2, trans_pars$`01`)
  timedeath <- rtrans(n, x1, x2, trans_pars$`02`)
  timesickdeath <- rtrans(n, x1, x2, trans_pars$`12`)
  censtime <- runif(n, 0, cens.par)
  # indicators
  D01 <- as.integer(timesick < pmin(timedeath, censtime))
  D02 <- as.integer(timedeath < pmin(timesick, censtime))
  D12 <- as.integer((D01 == 1) & ((timesick + timesickdeath) < censtime))
  # censored transition times
  obstimesick <- pmin(timesick, timedeath, censtime)
  obstimedeath <- ifelse(D02 == 1, timedeath,ifelse(D12 == 1, timesick + timesickdeath,censtime))
  # death indicator 
  id <- 1:n
  D2 <- D02 + D01 * D12
  # wide format
  simdat <- data.frame(id = id,obstimesick = obstimesick,D01 = D01,
    obstimedeath = obstimedeath,D2 = D2,x1 = x1,x2 = x2)
  # trnas matrix
  tmat <- transMat(x = list(c(2, 3), c(3), c()),names = c("Healthy", "Sick", "Dead"))
  mssimdat <- msprep(
    data = simdat,trans = tmat,
    time = c(NA, "obstimesick", "obstimedeath"),
    status = c(NA, "D01", "D2"),
    keep = c("x1", "x2")
  )
  mssimdat <- expand.covs(mssimdat,c("x1", "x2"),longnames = TRUE,append = TRUE)
  return(mssimdat)
}

eval_X_nl <- function( datas, ...) {
  
  mssimdat <- datas
  
  N <- length(unique(mssimdat$id))
  ind.test <- unique(mssimdat$id)
  
  Ni <- unlist(lapply(1:N, function(j, ...) {
    length(mssimdat$id[mssimdat$id == ind.test[j]])
  }))
  
  x1 <- matrix(NA, N, max(Ni))
  x2 <- matrix(NA, N, max(Ni))
  
  time.test <- matrix(NA, N, max(Ni))
  cen.test <- matrix(NA, N, max(Ni))
  trans.test <- matrix(NA, N, max(Ni))
  id.test <- matrix(NA, N, max(Ni))
  
  time.is.new <- mssimdat$time
  cens <- mssimdat$time
  time.is.new[mssimdat$status == 0] <- NA
  is.censored.new <- as.numeric(is.na(time.is.new))
  is.censored.test <- matrix(NA, N, max(Ni))
  
  count <- 1
  for (i in 1:N) {
    idx <- count:(Ni[i] + count - 1)
    
    x1[i, 1:Ni[i]] <- mssimdat[idx, "x1"]
    x2[i, 1:Ni[i]] <- mssimdat[idx, "x2"]
    
    time.test[i, 1:Ni[i]] <- time.is.new[idx]
    cen.test[i, 1:Ni[i]] <- cens[idx]
    trans.test[i, 1:Ni[i]] <- mssimdat[idx, "trans"]
    id.test[i, 1:Ni[i]] <- mssimdat[idx, "id"]
    is.censored.test[i, 1:Ni[i]] <- is.censored.new[idx]
    
    count <- count + Ni[i]
  }
  
  X <- model.matrix(~ x1[,1] + as.factor(x2[,1]))
  
  lista <- list(
    mssimdat = mssimdat,
    N = N,
    time.test = time.test,
    cen.test = cen.test,
    Ni = Ni,
    X = X,
    is.censored.test = is.censored.test,
    trans.test = trans.test,
    det = det(solve(t(X) %*% X)) != 0
  )
  
  return(lista)
}


# simulation structure ----------------------------------------------------
ns <- c(250,500,1000)
#* pars5: non linear ----------------------------------------------------
scenario_nonlinear <- list(`01` = list(type = "nonlinear", lambda = 1.2, k = 1.6, beta1 =  3, beta2 = -0.8),
                           `02` = list(type = "linear",    lambda = 0.6, k = 1.0, beta1 =  0.4, beta2 =  0.6),
                           `12` = list(type = "linear",    lambda = 1.0, k = 1.4, beta1 = -0.3, beta2 =  0.7) )

#* pars6: non ph -----------------------------------------------
scenario_nph <- list(`01` = list(type = "nph_logt",  lambda = 0.8, k = 1.2, beta1 =  0.8, beta2 = -0.3),
                     `02` = list(type = "linear",    lambda = 0.6, k = 1.0, beta1 =  0.4, beta2 =  0.6),
                     `12` = list(type = "nph_logt",  lambda = 1.0, k = 1.4, beta1 = -0.3, beta2 =  0.2) )

# * theta5 ----------------------------------------------------------------
# data: n1-par5 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data_nl(ns[1],trans_pars = scenario_nonlinear,cens.par = 2)
}
datas <- lapply(1:100,function(h) eval_X_nl(h,datas = data_trans[[h]]))

# data: n2-par5 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data_nl(ns[2],trans_pars = scenario_nonlinear,cens.par = 2)
}
datas <- lapply(1:100,function(h) eval_X_nl(h,datas = data_trans[[h]]))

# data: n3-par5 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data_nl(ns[3],trans_pars = scenario_nonlinear,cens.par = 2)
}
datas <- lapply(1:100,function(h) eval_X_nl(h,datas = data_trans[[h]]))

# * theta6 ----------------------------------------------------------------
# data: n1-par6 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data_nl(ns[1],trans_pars = scenario_nph,cens.par = 2)
}
datas <- lapply(1:100,function(h) eval_X_nl(h,datas = data_trans[[h]]))

# data: n2-par6 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data_nl(ns[2],trans_pars = scenario_nph,cens.par = 2)
}
datas <- lapply(1:100,function(h) eval_X_nl(h,datas = data_trans[[h]]))

# data: n3-par6 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data_nl(ns[3],trans_pars = scenario_nph,cens.par = 2)
}
datas <- lapply(1:100,function(h) eval_X_nl(h,datas = data_trans[[h]]))






# ---------------------------------------------------------
# 3. Fit multistate clock-reset Cox model
# ---------------------------------------------------------
fit_ms_cox <- function(mssimdat) {
  coxph(
    Surv(time, status) ~ strata(trans) + x1.1 + x1.2 + x1.3 + x2.1 + x2.2 + x2.3,
    data = mssimdat
  )
}

# ---------------------------------------------------------
# 4. Extract baseline cumulative hazard by transition
# ---------------------------------------------------------
get_cox_hazards <- function(cox_fit, mssimdat) {
  Haz.c <- basehaz(cox_fit, centered = FALSE)
  all_transitions <- sort(unique(mssimdat$trans))
  
  Haz_strata <- lapply(all_transitions, function(tr) {
    Haz.c[Haz.c$strata == paste0("trans=",tr), ]
  })
  
  names(Haz_strata) <- as.character(all_transitions)
  Haz_strata
}

# ---------------------------------------------------------
# 5. Predict Cox survival at Cox jump times for one row
# ---------------------------------------------------------
predict_cox_surv_at_jump_times <- function(cox_fit, haz_list, nd_row, tr) {
  haz <- haz_list[[as.character(tr)]]
  
  t_cox <- haz$time
  H0 <- haz$hazard
  
  lin_pred <- sum(coef(cox_fit) * as.numeric(nd_row[names(coef(cox_fit))]))
  surv_cox <- exp(-H0 * exp(lin_pred))
  
  list(time = t_cox, surv = surv_cox)
}

# ---------------------------------------------------------
# 6. Evaluate one subject under one scenario
# ---------------------------------------------------------
evaluate_subject_cox_general <- function(h, mssimdat, cox_fit, haz_list, scenario) {
  
  nd <- mssimdat[mssimdat$id == h, ]
  x1_val <- unique(nd$x1)
  x2_val <- unique(nd$x2)
  trans_numbs <- as.numeric(nd$trans)
  
  difs_subject <- matrix(0, nrow = length(trans_numbs), ncol = 2)
  colnames(difs_subject) <- c("L1", "L2")
  
  trans_map <- c("01", "02", "12")
  
  for (i in seq_along(trans_numbs)) {
    tr <- trans_numbs[i]
    tr_name <- trans_map[tr]
    
    pred <- predict_cox_surv_at_jump_times(
      cox_fit = cox_fit,
      haz_list = haz_list,
      nd_row = nd[i, ],
      tr = tr
    )
    
    t_cox <- pred$time
    surv_cox <- pred$surv
    
    real <- true_surv(
      t = t_cox,
      pars = scenario[[tr_name]],
      x1 = x1_val,
      x2 = x2_val
    )
    
    difs_subject[i, ] <- c(
      sum(abs(real - surv_cox)),
      sum((real - surv_cox)^2)
    )
  }
  
  colSums(difs_subject)
}

# ---------------------------------------------------------
# 7. Evaluate one full dataset
# ---------------------------------------------------------
evaluate_dataset_cox_general <- function(mssimdat, scenario) {
  cox_fit <- fit_ms_cox(mssimdat)
  haz_list <- get_cox_hazards(cox_fit, mssimdat)
  
  subject_ids <- unique(mssimdat$id)
  
  out <- t(sapply(subject_ids, function(h) {
    evaluate_subject_cox_general(
      h = h,
      mssimdat = mssimdat,
      cox_fit = cox_fit,
      haz_list = haz_list,
      scenario = scenario
    )
  }))
  
  out <- as.data.frame(out)
  out$id <- subject_ids
  
  list(
    cox_fit = cox_fit,
    hazards = haz_list,
    subject_results = out,
    mean_L1 = mean(out$L1),
    mean_L2 = mean(out$L2)
  )
}


# ---------------------------------------------------------
# 9. EXAMPLES OF USE
# ---------------------------------------------------------
# IMPORTANT:
# You must already have mssimdat for theta5/theta6.
# For example:
n1.mssimdat_theta5 <- lapply(1:100,function(h) {eval_X_nl( sim_data_nl(n = 250, trans_pars = scenario_nonlinear, cens.par = 2) )})
n2.mssimdat_theta5 <- lapply(1:100,function(h) {eval_X_nl( sim_data_nl(n = 500, trans_pars = scenario_nonlinear, cens.par = 2) )})
n3.mssimdat_theta5 <- lapply(1:100,function(h) {eval_X_nl( sim_data_nl(n = 1000, trans_pars = scenario_nonlinear, cens.par = 2) )})

# Example for theta5:
n1.res_theta5 <- lapply(1:100,function(h){ evaluate_dataset_cox_general(
  mssimdat = n1.mssimdat_theta5[[h]]$mssimdat,
  scenario = scenario_nonlinear
)})

n2.res_theta5 <- lapply(1:100,function(h){ evaluate_dataset_cox_general(
  mssimdat = n2.mssimdat_theta5[[h]]$mssimdat,
  scenario = scenario_nonlinear
)})

n3.res_theta5 <- lapply(1:100,function(h){ evaluate_dataset_cox_general(
  mssimdat = n3.mssimdat_theta5[[h]]$mssimdat,
  scenario = scenario_nonlinear
)})

# You must already have mssimdat for theta6/theta6.
# For example:
n1.mssimdat_theta6 <- lapply(1:100,function(h) {eval_X_nl( sim_data_nl(n = 250, trans_pars = scenario_nph, cens.par = 2) )})
n2.mssimdat_theta6 <- lapply(1:100,function(h) {eval_X_nl( sim_data_nl(n = 500, trans_pars = scenario_nph, cens.par = 2) )})
n3.mssimdat_theta6 <- lapply(1:100,function(h) {eval_X_nl( sim_data_nl(n = 1000, trans_pars = scenario_nph, cens.par = 2) )})

# Example for theta6:
n1.res_theta6 <- lapply(1:100,function(h){ evaluate_dataset_cox_general(
  mssimdat = n1.mssimdat_theta6[[h]]$mssimdat,
  scenario = scenario_nph
)})

n2.res_theta6 <- lapply(1:100,function(h){ evaluate_dataset_cox_general(
  mssimdat = n2.mssimdat_theta6[[h]]$mssimdat,
  scenario = scenario_nph
)})

n3.res_theta6 <- lapply(1:100,function(h){ evaluate_dataset_cox_general(
  mssimdat = n3.mssimdat_theta6[[h]]$mssimdat,
  scenario = scenario_nph
)})


mean(unlist(lapply(1:100, function(h) n1.res_theta5[[h]]$mean_L1)))
mean(unlist(lapply(1:100, function(h) n2.res_theta5[[h]]$mean_L1)))
mean(unlist(lapply(1:100, function(h) n3.res_theta5[[h]]$mean_L1)))

mean(unlist(lapply(1:100, function(h) n1.res_theta5[[h]]$mean_L2)))
mean(unlist(lapply(1:100, function(h) n2.res_theta5[[h]]$mean_L2)))
mean(unlist(lapply(1:100, function(h) n3.res_theta5[[h]]$mean_L2)))


mean(unlist(lapply(1:100, function(h) n1.res_theta6[[h]]$mean_L1)))
mean(unlist(lapply(1:100, function(h) n2.res_theta6[[h]]$mean_L1)))
mean(unlist(lapply(1:100, function(h) n3.res_theta6[[h]]$mean_L1)))

mean(unlist(lapply(1:100, function(h) n1.res_theta6[[h]]$mean_L2)))
mean(unlist(lapply(1:100, function(h) n2.res_theta6[[h]]$mean_L2)))
mean(unlist(lapply(1:100, function(h) n3.res_theta6[[h]]$mean_L2)))



mssimdat_theta6 <- sim_data_nl(n = 500, trans_pars = scenario_nph, cens.par = 2)
