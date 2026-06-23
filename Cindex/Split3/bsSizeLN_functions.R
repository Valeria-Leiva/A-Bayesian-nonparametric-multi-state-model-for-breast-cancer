# Evaluation predictive distributions
# by transitions - hospital - covars

# rm(list=ls())


# set working directory ---------------------------------------------------
setwd("..../Split3")


# libraries ---------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse,survival,ggplot2,splines)

library("stats")

 # BNP fit -----------------------------------------------------------------


options("max.print"=1000000)
mem.maxVSize(32000000000)

load("Res_bsSizeLNO_BNP_1_zellnerTrain.RData")
betah2 <- out.DDP2$BUGSoutput$sims.list$betah
# dim(betah2) #Nmc p Tr Hdp
bk2 <- out.DDP2$BUGSoutput$sims.list$bk
# dim(bk2) #Nmc hosp
wh2 <- out.DDP2$BUGSoutput$sims.list$BetaProb
# dim(wh2) #Nmc Hdp Tr
sig2 <- out.DDP2$BUGSoutput$sims.list$sigma2
# dim(sig2) #Nmc  Tr
M2<- out.DDP2$BUGSoutput$sims.list$M
# dim(M2) #Nmc 1
nstar <- out.DDP2$BUGSoutput$sims.list$m
# dim(nstar) #Nmc n Tr


Nmc <- dim(betah2)[1]
Hdp <- dim(betah2)[4]
tra <- dim(betah2)[3]
p <- dim(betah2)[2]
# nper <- dim(nstar)[2]

pred_dist <- function(y,xx,tr,hp){
  idh = hp
  Fmcdp <- NULL
  for(imc in 1:Nmc){
    Fmcdp[imc]=1-rowSums( do.call('cbind',lapply(1:Hdp,function(jdp) wh2[imc,jdp,tr]*plnorm(y,meanlog=as.numeric(betah2[imc,,tr,jdp]%*%xx)+bk2[imc,idh],
                                                                                            sdlog = sqrt(sig2[imc,tr])))))
  }
  return(Fmcdp)
}

eval_p=function(y,xx,tr,hp){
  eval.pred = lapply(1:length(y),function(h) pred_dist(y[h],xx,tr,hp))
  eval.pred.mean = unlist(lapply(eval.pred,mean))
  eval.pred.q =  do.call('rbind',lapply(eval.pred,function(h) quantile(h,c(0.025,0.975))))
  return(cbind(
    pred_mean=eval.pred.mean,
    pred_q=eval.pred.q
  ))
}

eval_p_mean=function(y,xx,tr,hp){
  eval.pred = lapply(1:length(y),function(h) pred_dist(y[h],xx,tr,hp))
  eval.pred.mean = unlist(lapply(eval.pred,mean))
  # eval.pred.q =  do.call('rbind',lapply(eval.pred,function(h) quantile(h,c(0.025,0.975))))
  return(cbind(
    pred_mean=eval.pred.mean
    # ,pred_q=eval.pred.q
  ))
}

eval_p_cindex=function(y.x,xx.x,tr.x,hp.x,
                       y.y,xx.y,tr.y,hp.y){
  eval.pred.x = lapply(1:length(y.x),function(h) pred_dist(y.x[h],xx.x,tr.x,hp.x))
  eval.pred.y = lapply(1:length(y.y),function(h) pred_dist(y.y[h],xx.y,tr.y,hp.y))
  
  return(mean(eval.pred.x[[1]]<eval.pred.y[[1]]))
}


