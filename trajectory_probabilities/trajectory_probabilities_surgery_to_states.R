## post-surgery

rm(list = ls())
setwd("...../H25")


library(survival)
library(ggplot2)
library(parallel)

load('h25.Rdata')

Nmc <- dim(betah2)[1]
Hdp <- dim(betah2)[4]
tra <- dim(betah2)[3]
p <- dim(betah2)[2]

# Predictive ---------------------------------------------------------------
#hosp=hosp.test[,1]
pred_dist <- function(y,xx,tr,hp){

  idh = hp
  #nper = sum(idh)
  #Fmcdp <- array(NA,dim=c(Nmc,nper))
  Fmcdp <- NULL
  for(imc in 1:Nmc){
    Fmcdp[imc]=1-rowSums( do.call('cbind',lapply(1:Hdp,function(jdp) wh2[imc,jdp,tr]*plnorm(y,meanlog=as.numeric(betah2[imc,,tr,jdp]%*%xx)+bk2[imc,idh],
                                                                                             sdlog = sqrt(sig2[imc,tr])))))
  }
  return(Fmcdp)
}
pred_hazard <- function(y,xx,tr,hp){

  idh = hp
  #nper = sum(idh)
  #Fmcdp <- array(NA,dim=c(Nmc,nper))
  Fmcdp <- NULL
  fmcdp <- NULL


  for(imc in 1:Nmc){
    Fmcdp[imc]=1-rowSums( do.call('cbind',lapply(1:Hdp,function(jdp) wh2[imc,jdp,tr]*plnorm(y,meanlog=as.numeric(betah2[imc,,tr,jdp]%*%xx)+bk2[imc,idh],
                                                                                            sdlog = sqrt(sig2[imc,tr])))))
    fmcdp[imc]=rowSums( do.call('cbind',lapply(1:Hdp,function(jdp) wh2[imc,jdp,tr]*dlnorm(y,meanlog=as.numeric(betah2[imc,,tr,jdp]%*%xx)+bk2[imc,idh],
                                                                                            sdlog = sqrt(sig2[imc,tr])))))

  }
  hmcdp <- fmcdp/Fmcdp
  return(hmcdp)
}
# Evaluation --------------------------------------------------------------

eval_p=function(y,xx,tr,hp){
 # eval.post = lapply(1:length(y),function(h) post_dist(y[h],xx,tr,hp))
 # eval.post.mean = unlist(lapply(eval.post,mean))
 # eval.post.q =  do.call('rbind',lapply(eval.post,function(h) quantile(h,c(0.025,0.975))))

  eval.pred = lapply(1:length(y),function(h) pred_dist(y[h],xx,tr,hp))
  eval.pred.mean = unlist(lapply(eval.pred,mean))
  eval.pred.q =  do.call('rbind',lapply(eval.pred,function(h) quantile(h,c(0.025,0.975))))
  return(list(
    #post_mean=eval.post.mean,
   # post_q=eval.post.q,
    pred_mean=eval.pred.mean,pred_q=eval.pred.q
  ))
}


eval_h=function(y,xx,tr,hp){

  eval.pred = lapply(1:length(y),function(h) pred_hazard(y[h],xx,tr,hp))
  eval.pred.mean = unlist(lapply(eval.pred,mean))
  eval.pred.q =  do.call('rbind',lapply(eval.pred,function(h) quantile(h,c(0.025,0.975))))
  return(list(
    pred_mean=eval.pred.mean,pred_q=eval.pred.q
  ))
}

## choose other
new_subject <- c(1,57,0,0.3657755,0.081596068,0.0060673999,1,0,1) ## ER(+)
#new_subject <- c(1,52,4,0.39648688,0.1078444324,9.777895e-03,0,1,0) ## ER(-)
timepoints <- seq(from=0, to=20, by=0.25)
#hp = 1


# Distant Relapse
S2 <- function(u,...) {
  # t : start
  #out1 <- -log(eval_p(u,new_subject,tr=8,hp=3)$pred_mean)- (-log(eval_p(t,new_subject,tr=8,hp=5)$pred_mean))
  #out2 <- -log(eval_p(u,new_subject,tr=9,hp=3)$pred_mean) - (-log(eval_p(t,new_subject,tr=9,hp=5)$pred_mean))
  #aux <- exp(-(out1+out2))
  out1 = eval_p(u,new_subject,tr=8,hp=hp)$pred_mean
  out2 = eval_p(u,new_subject,tr=9,hp=hp)$pred_mean
  aux = out1*out2
  return(aux)
}
## Local Relapse
S1 <- function(u,...) {
  # t : start
  #out1 <- -log(eval_p(u,new_subject,tr=8,hp=3)$pred_mean)- (-log(eval_p(t,new_subject,tr=8,hp=5)$pred_mean))
  #out2 <- -log(eval_p(u,new_subject,tr=9,hp=3)$pred_mean) - (-log(eval_p(t,new_subject,tr=9,hp=5)$pred_mean))
  #aux <- exp(-(out1+out2))
  out1 = eval_p(u,new_subject,tr=5,hp=hp)$pred_mean
  out2 = eval_p(u,new_subject,tr=6,hp=hp)$pred_mean
  out3 = eval_p(u,new_subject,tr=7,hp=hp)$pred_mean
  aux = out1*out2*out3
  return(aux)
}





S3 <- function(u,...) {
  # t : start
  #out1 <- -log(eval_p(u,new_subject,tr=8,hp=3)$pred_mean)- (-log(eval_p(t,new_subject,tr=8,hp=5)$pred_mean))
  #out2 <- -log(eval_p(u,new_subject,tr=9,hp=3)$pred_mean) - (-log(eval_p(t,new_subject,tr=9,hp=5)$pred_mean))
  #aux <- exp(-(out1+out2))
  out1 = eval_p(u,new_subject,tr=1,hp=hp)$pred_mean
  out2 = eval_p(u,new_subject,tr=2,hp=hp)$pred_mean
  out3 = eval_p(u,new_subject,tr=3,hp=hp)$pred_mean
  out4 = eval_p(u,new_subject,tr=4,hp=hp)$pred_mean
  aux = out1*out2*out3*out4
  return(aux)
}



## SO
target_pso <- function(x,...){
  nintervals = (x - 0)/20
  intervals = seq(0, x, by = nintervals)
  interval.widths = intervals[2:length(intervals)] - intervals[1:(length(intervals)-1)]
  midpoints = 0.5*(intervals[2:length(intervals)] + intervals[1:(length(intervals)-1)])
  aux = sapply(1:length(midpoints),function(j,...){
  	(eval_h(midpoints[j],new_subject,tr=4,hp=hp)$pred_mean*S3(midpoints[j]))*interval.widths[j]
  })
  return(sum(aux))
}
## SC
target_psc <- function(x,...){
  nintervals = (x - 0)/20
  intervals = seq(0, x, by = nintervals)
  interval.widths = intervals[2:length(intervals)] - intervals[1:(length(intervals)-1)]
  midpoints = 0.5*(intervals[2:length(intervals)] + intervals[1:(length(intervals)-1)])
  aux = sapply(1:length(midpoints), function(j,...){
  	(eval_h(midpoints[j],new_subject,tr=3,hp=hp)$pred_mean*S3(midpoints[j]))*interval.widths[j]
  })
  return(sum(aux))
}

## SDC
target_psdc <- function(x,...){
  n_intervals <- 20
  step <- x / n_intervals	
  # Create intervals and midpoints
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)
  
  # Precompute expensive function values
  eval_tr2 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 2, hp = hp)$pred_mean)
  eval_tr8 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 8, hp = hp)$pred_mean)
  S2_vals <- sapply(midpoints, S2)
  S3_vals <- sapply(midpoints, S3)	
  eval_tr8_rev <- rev(eval_tr8)
  S2_vals_rev <- rev(S2_vals)
  widths_rev <- rev(widths)
  
  # Compute outer matrix
  outer_matrix <- outer(
    1:n_mid, 1:n_mid,
    Vectorize(function(j, k) {
      eval_tr2[j] * S3_vals[j] * widths[j] *
        eval_tr8_rev[k] * S2_vals_rev[k] * widths_rev[k]
    })
  )

  return(sum(outer_matrix))

}
## SDO
target_psdo <- function(x,...){
  n_intervals <- 20
  step <- x / n_intervals	
  # Create intervals and midpoints
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)

    # Precompute expensive function values
  eval_tr2 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 2, hp = hp)$pred_mean)
  eval_tr9 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 9, hp = hp)$pred_mean)
  S2_vals <- sapply(midpoints, S2)
  S3_vals <- sapply(midpoints, S3)	
  eval_tr9_rev <- rev(eval_tr9)
  S2_vals_rev <- rev(S2_vals)
  widths_rev <- rev(widths)

  # Compute outer matrix
  outer_matrix <- outer(
    1:n_mid, 1:n_mid,
    Vectorize(function(j, k) {
      eval_tr2[j] * S3_vals[j] * widths[j] *
        eval_tr9_rev[k] * S2_vals_rev[k] * widths_rev[k]
    })
  )

  return(sum(outer_matrix))

}



## SLO
target_pslo <- function(x,...){
  n_intervals <- 20
  step <- x / n_intervals	
  # Create intervals and midpoints
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)
  
      # Precompute expensive function values
  eval_tr1 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 1, hp = hp)$pred_mean)
  eval_tr7 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 7, hp = hp)$pred_mean)
  S1_vals <- sapply(midpoints, S1)
  S3_vals <- sapply(midpoints, S3)	
  eval_tr7_rev <- rev(eval_tr7)
  S1_vals_rev <- rev(S1_vals)
  widths_rev <- rev(widths)
  
  # Compute outer matrix
  outer_matrix <- outer(
    1:n_mid, 1:n_mid,
    Vectorize(function(j, k) {
      eval_tr1[j] * S3_vals[j] * widths[j] *
        eval_tr7_rev[k] * S1_vals_rev[k] * widths_rev[k]
    })
  )

  return(sum(outer_matrix))

}


## SLC
target_pslc <- function(x,...){
n_intervals <- 20
  step <- x / n_intervals	
  # Create intervals and midpoints
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)
  
   # Precompute expensive function values
  eval_tr1 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 1, hp = hp)$pred_mean)
  eval_tr6 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 6, hp = hp)$pred_mean)
  S1_vals <- sapply(midpoints, S1)
  S3_vals <- sapply(midpoints, S3)	
  eval_tr6_rev <- rev(eval_tr6)
  S1_vals_rev <- rev(S1_vals)
  widths_rev <- rev(widths)


 # Compute outer matrix
  outer_matrix <- outer(
    1:n_mid, 1:n_mid,
    Vectorize(function(j, k) {
      eval_tr1[j] * S3_vals[j] * widths[j] *
        eval_tr6_rev[k] * S1_vals_rev[k] * widths_rev[k]
    })
  )

  return(sum(outer_matrix))
}


## SLDO


target_psldo <- function(x, ...) {
  n_intervals <- 20
  step <- x / n_intervals

  # Compute intervals and midpoints
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)

  # Precompute values outside the loop
  eval_tr1 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 1, hp = hp)$pred_mean)
  eval_tr5 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 5, hp = hp)$pred_mean)
  eval_tr9 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 9, hp = hp)$pred_mean)
  S1_vals <- sapply(midpoints, S1)
  S2_vals <- sapply(midpoints, S2)
  S3_vals <- sapply(midpoints, S3)

  # Reverse vectors for k-indexed operations
  eval_tr5_rev <- rev(eval_tr5)
  eval_tr9_rev <- rev(eval_tr9)
  S1_rev <- rev(S1_vals)
  S2_rev <- rev(S2_vals)
  widths_rev <- rev(widths)

  # Precompute auxdo and auxd (length = n_mid)
  auxdo <- eval_tr9_rev * S2_rev * widths_rev
  auxd <- eval_tr5_rev * auxdo * S1_rev

  # Compute outer product
  outer_result <- outer(
    1:n_mid, 1:n_mid,
    Vectorize(function(j, k) {
      eval_tr1[j] * auxdo[k] * auxd[k] * S3_vals[j] * widths[j]
    })
  )

  return(sum(outer_result))
}


## SLDC
target_psldo <- function(x, ...) {
  n_intervals <- 20
  step <- x / n_intervals

  # Compute intervals and midpoints
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)

  # Precompute values outside the loop
  eval_tr1 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 1, hp = hp)$pred_mean)
  eval_tr5 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 5, hp = hp)$pred_mean)
  eval_tr8 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 8, hp = hp)$pred_mean)
  S1_vals <- sapply(midpoints, S1)
  S2_vals <- sapply(midpoints, S2)
  S3_vals <- sapply(midpoints, S3)

  # Reverse vectors for k-indexed operations
  eval_tr5_rev <- rev(eval_tr5)
  eval_tr8_rev <- rev(eval_tr8)
  S1_rev <- rev(S1_vals)
  S2_rev <- rev(S2_vals)
  widths_rev <- rev(widths)

  # Precompute auxdo and auxd (length = n_mid)
  auxdo <- eval_tr8_rev * S2_rev * widths_rev
  auxd <- eval_tr5_rev * auxdo * S1_rev

  # Compute outer product
  outer_result <- outer(
    1:n_mid, 1:n_mid,
    Vectorize(function(j, k) {
      eval_tr1[j] * auxdo[k] * auxd[k] * S3_vals[j] * widths[j]
    })
  )

  return(sum(outer_result))
}


## SLD

target_psld <- function(x, ...) {
  n_intervals <- 20
  step <- x / n_intervals
  
  # Construct intervals, midpoints, and widths
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)

  # Precompute all function evaluations
  eval_tr1 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 1, hp = hp)$pred_mean)
  eval_tr5 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 5, hp = hp)$pred_mean)
  eval_tr8 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 8, hp = hp)$pred_mean)
  eval_tr9 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 9, hp = hp)$pred_mean)
  
  S1_vals <- sapply(midpoints, S1)
  S2_vals <- sapply(midpoints, S2)
  S3_vals <- sapply(midpoints, S3)

  # Reverse for k-indexed values
  eval_tr5_rev <- rev(eval_tr5)
  eval_tr8_rev <- rev(eval_tr8)
  eval_tr9_rev <- rev(eval_tr9)
  S1_rev <- rev(S1_vals)
  S2_rev <- rev(S2_vals)
  widths_rev <- rev(widths)

  # Precompute terms for each component
  auxdc <- eval_tr8_rev * S2_rev * widths_rev
  auxdc1 <- eval_tr5_rev * auxdc * S1_rev * widths_rev

  auxdo <- eval_tr9_rev * S2_rev * widths_rev
  auxdo1 <- eval_tr5_rev * auxdo * S1_rev * widths_rev

  auxdl <- eval_tr5_rev * S1_rev * widths_rev

  # Precompute j-side values
  outer_factor <- eval_tr1 * S3_vals * widths

  # Compute all three matrices using outer product
  aux1 <- outer(outer_factor, auxdc1)
  aux2 <- outer(outer_factor, auxdo1)
  aux3 <- outer(outer_factor, auxdl)

  # Final result
  return(sum(aux3) - sum(aux1) - sum(aux2))
}


## SD
target_psd <- function(x, ...) {
  n_intervals <- 20
  step <- x / n_intervals

  # Midpoints and interval widths
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, n_intervals)
  n_mid <- length(midpoints)

  # Precompute all eval_h and S functions
  eval_tr2 <- sapply(midpoints, function(t) eval_h(t, new_subject, tr = 2, hp = hp)$pred_mean)
  eval_tr8 <- sapply(rev(midpoints), function(t) eval_h(t, new_subject, tr = 8, hp = hp)$pred_mean)
  eval_tr9 <- sapply(rev(midpoints), function(t) eval_h(t, new_subject, tr = 9, hp = hp)$pred_mean)

  S2_rev <- sapply(rev(midpoints), S2)
  S3_vals <- sapply(midpoints, S3)
  widths_rev <- rev(widths)

  # Precompute inner integrals
  auxdc <- eval_tr8 * S2_rev * widths_rev
  auxdo <- eval_tr9 * S2_rev * widths_rev

  # Compute outer product terms
  outer_factor <- eval_tr2 * S3_vals * widths
  aux1 <- outer(outer_factor, auxdc)
  aux2 <- outer(outer_factor, auxdo)
  aux0 <- outer_factor

  # Final computation
  return(sum(aux0) - sum(aux1) - sum(aux2))
}


## SL


target_psl <- function(x, ...) {
  n_intervals <- 20
  step <- x / n_intervals

  # Set up midpoints and interval widths
  intervals <- seq(0, x, by = step)
  midpoints <- 0.5 * (head(intervals, -1) + tail(intervals, -1))
  widths <- rep(step, length(midpoints))
  rev_midpoints <- rev(midpoints)
  rev_widths <- rev(widths)

  # Precompute evaluations (memoization)
  evals <- function(tr, times) {
    sapply(times, function(t) eval_h(t, new_subject, tr = tr, hp = hp)$pred_mean)
  }

  S1_rev <- sapply(rev_midpoints, S1)
  S2_rev <- sapply(rev_midpoints, S2)
  S3_vals <- sapply(midpoints, S3)

  # Precompute outer eval_h
  eval_1 <- evals(1, midpoints)
  eval_2 <- evals(2, midpoints)

  # Precompute inner eval_h
  eval_5_rev <- evals(5, rev_midpoints)
  eval_6_rev <- evals(6, rev_midpoints)
  eval_7_rev <- evals(7, rev_midpoints)
  eval_8_rev <- evals(8, rev_midpoints)
  eval_9_rev <- evals(9, rev_midpoints)

  # Precompute components
  auxdc  <- eval_8_rev * S2_rev * rev_widths
  auxdo  <- eval_9_rev * S2_rev * rev_widths
  auxld  <- eval_5_rev * S1_rev * rev_widths
  auxldc <- eval_5_rev * auxdc  * S1_rev * rev_widths
  auxldo <- eval_5_rev * auxdo  * S1_rev * rev_widths
  auxlc  <- eval_6_rev * S1_rev * rev_widths
  auxlo  <- eval_7_rev * S1_rev * rev_widths

  # Prepare outer weights
  outer_weight_1 <- eval_1 * S3_vals * widths
  outer_weight_2 <- eval_2 * S3_vals * widths

  # Use outer product for integration
  aux1 <- outer(outer_weight_1, eval_5_rev * auxdc * S1_rev * rev_widths)
  aux2 <- outer(outer_weight_1, eval_5_rev * auxdo * S1_rev * rev_widths)
  aux3 <- outer(outer_weight_1, auxld - auxldc - auxldo)
  aux4 <- outer(outer_weight_1, auxlc)
  aux5 <- outer(outer_weight_2, auxlo)

  aux0 <- outer_weight_1

  return(sum(aux0) - sum(aux1) - sum(aux2) - sum(aux3) - sum(aux4) - sum(aux5))
}



evaluate_hp <- function(hp_value){
	assign('hp', hp_value, envir = .GlobalEnv)
	n.t <- length(timepoints)
	prob_pso0 <- numeric(n.t)
	prob_psc0 <- numeric(n.t)
	prob_psdc0 <- numeric(n.t)
	prob_psdo0 <- numeric(n.t)
	prob_pslo0 <- numeric(n.t)
	prob_pslc0 <- numeric(n.t)
	prob_psldo0 <- numeric(n.t)
	prob_psldc0 <- numeric(n.t)
	prob_psld0 <- numeric(n.t)
	prob_psd0 <- numeric(n.t)
	prob_psl0 <- numeric(n.t)

	prob_pso0[1] <- 0
	prob_psc0[1] <- 0
	prob_psdc0[1] <- 0
	prob_psdo0[1] <- 0
	prob_pslo0[1] <- 0
	prob_pslc0[1] <- 0
	prob_psldo0[1] <- 0
	prob_psldc0[1] <- 0
	prob_psld0[1] <- 0
	prob_psd0[1] <- 0
	prob_psl0[1] <- 0
	
	indices <- 2:n.t	
	
	prob_pso0[indices] <- sapply(timepoints[indices], target_pso)
	prob_psc0[indices] <- sapply(timepoints[indices], target_psc)
	prob_psdc0[indices] <- sapply(timepoints[indices], target_psdc)
	prob_psdo0[indices] <- sapply(timepoints[indices],target_psdo)
	prob_pslo0[indices] <- sapply(timepoints[indices], target_pslo)
	prob_pslc0[indices] <- sapply(timepoints[indices], target_pslc)
	prob_psldo0[indices] <- sapply(timepoints[indices], target_psldo)
	prob_psldc0[indices] <- sapply(timepoints[indices], target_psldc)
	prob_psld0[indices] <- sapply(timepoints[indices], target_psld)
	prob_psd0[indices] <- sapply(timepoints[indices], target_psd)
	prob_psl0[indices] <- sapply(timepoints[indices], target_psl)
	
  data.frame('hp'= hp_value, 'time' = timepoints,
  'prob_pso' = prob_pso0, 'prob_psc'= prob_psc0, 'prob_psdc'= prob_psdc0, 
'prob_psdo'= prob_psdo0,'prob_pslo'= prob_pslo0,
 'prob_pslc'= prob_pslc0, 'prob_psldo'= prob_psldo0,
 'prob_psldc' = prob_psldc0, 'prob_psld'= prob_psld0,'prob_psd'= prob_psd0,
 'prob_psl' = prob_psl0)
 
}

results_list <- mclapply(1:5, evaluate_hp, mc.cores = 5)
results_df <- do.call(rbind, results_list)


save(results_df, file = 'Probs_surgery_posh.Rdata')
