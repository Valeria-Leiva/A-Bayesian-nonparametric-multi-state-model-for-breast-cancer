#simulation study: functions

# libraries ---------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel,future,furrr,
               tidyverse,survival,splines,mstate,
               rjags,R2jags)


# functions ---------------------------------------------------------------
#simulation data
sim_data <- function(n,beta,rate,cens.par,...){
  x1 <- runif(n)
  x2 <- rbinom(n,size = 1,prob = 0.5)
  u <- runif(n)
  timesick <- ( -log(1-u)/( rate[1]* exp(beta[1]*x1 + beta[2]*x2) ) )^(1/ rate[2]) 
  u <- runif(n)
  timedeath <- ( -log(1-u)/ ( rate[3]* exp(beta[3]*x1+ beta[4]*x2) ) )^(1/ rate[4])
  u <- runif(n)
  timesickdeath <- ( -log(1-u)/ ( rate[5]* exp(beta[5]*x1 + beta[6]*x2) ) )^(1/ rate[6]) 
  censtime <- runif(n, 0, cens.par)
  # indicators
  D01 <- 1*(timesick < pmin(timedeath,censtime))
  D02 <- 1*(timedeath < pmin(timesick,censtime))
  D12 <- 1*(D01==1)*((timesick+timesickdeath)< censtime)
  # censored transition times
  obstimesick <- timesick*D01 + pmin(timedeath,censtime)*(1-D01)
  obstimedeath <- timedeath*D02 + censtime*D01*(1-D12) + (timesick+timesickdeath)*D12 + censtime*(1-D01-D02)
  #  death indicator
  id <- 1:n
  D2 <- D02+D01*D12
  #wide format 
  simdat <- data.frame(cbind(id,obstimesick,D01,obstimedeath,D2, x1, x2))
  # trnas matrix
  tmat <- transMat(x = list(c(2,3),c(3),c()),names = c("Healthy","Sick","Dead"))
  mssimdat <- msprep(data = simdat, trans=tmat,
                     time = c(NA,"obstimesick","obstimedeath"),
                     status = c(NA,"D01","D2"),
                     keep = c("x1","x2"))
  mssimdat<-expand.covs(mssimdat,c("x1","x2"),longnames = TRUE, append = TRUE)
  return(mssimdat)
}

#eval design matrix
eval_X <- function(h,datas,...){
  
  mssimdat <- datas
  
  N <- length(unique(mssimdat$id))
  ind.test <- unique(mssimdat$id)
  Ni <- unlist(lapply(1:N, function(j,...){
    length(mssimdat$id[mssimdat$id==ind.test[j]])
  }))
  
  x1 <- matrix(NA,N,max(Ni))
  x2 <- matrix(NA,N,max(Ni))
  x3 <- matrix(NA,N,max(Ni))
  
  time.test <- matrix(NA,N,max(Ni))
  cen.test <- matrix(NA,N,max(Ni))
  trans.test <- matrix(NA,N,max(Ni))
  id.test <- matrix(NA,N,max(Ni))
  
  
  ## solo para crear is.censored
  time.is.new<- mssimdat$time
  cens <- mssimdat$time
  time.is.new[mssimdat$status==0] <- NA
  is.censored.new <- as.numeric(is.na(time.is.new))
  is.censored.test <- matrix(NA,N,max(Ni))
  
  count <- 1
  for(i in 1:N){
    x1[i,1:Ni[i]] <- mssimdat[count:(Ni[i]+count-1),"x1"]
    x2[i,1:Ni[i]] <- mssimdat[count:(Ni[i]+count-1),"x2"]
    x3[i,1:Ni[i]] <- mssimdat[count:(Ni[i]+count-1),"x3"]
    time.test[i,1:Ni[i]] <- time.is.new[count:(Ni[i]+count-1)]
    cen.test[i,1:Ni[i]] <- cens[count:(Ni[i]+count-1)]
    trans.test[i,1:Ni[i]] <- mssimdat[count:(Ni[i]+count-1),"trans"]
    id.test[i,1:Ni[i]] <- mssimdat[count:(Ni[i]+count-1),"id"]
    is.censored.test[i,1:Ni[i]] <- is.censored.new[count:(Ni[i]+count-1)]
    count <- count + Ni[i]
  }
  
  X <- model.matrix(~x1[,1]+as.factor(x2[,1])+x3[,1])
  
  lista <- list(mssimdat = mssimdat,
                h=h, 
                N = N,
                time.test=time.test,
                cen.test = cen.test, 
                Ni = Ni,
                X = X, 
                is.censored.test = is.censored.test,
                trans.test = trans.test,
                det = det(solve(t(X)%*%X)) != 0)
  return(lista)
}
#evaluation real survival function
S_real <- function(tt,trans,betas,rates,x,...){
  
  i <- c(1,3,5)
  ii <- i[trans]
  
  eval <- exp(-(tt^rates[ii+1])*rates[ii]*exp(betas[ii]*x[1]+betas[ii+1]*x[2]))
  
  return(eval)
}
#evaluation predictive distribution
S_pred <- function(y,xx,tr,..){
  Fmcdp <- NULL
  for(imc in 1:Nmc){
    Fmcdp[imc]=1-rowSums( do.call('cbind',lapply(1:Hdp,function(jdp) wh2[imc,jdp,tr]*plnorm(y,meanlog=sum(betah2[imc,,tr,jdp]*xx),
                                                                                            sdlog = sqrt(sig2[imc,tr])))))
  }
  return(Fmcdp)
}
#summary predictive distribution
S_pred_summary <- function(y,xx,tr){
  eval.pred = lapply(1:length(y),function(h) pred_dist(y[h],xx,tr))
  eval.pred.mean = unlist(lapply(eval.pred,mean))
  eval.pred.q =  do.call('rbind',lapply(eval.pred,function(h) quantile(h,c(0.025,0.975))))
  return(list(
    pred_mean=eval.pred.mean,pred_q=eval.pred.q
  ))
}

S_pred_mean=function(y,xx,tr){
  eval.pred = lapply(1:length(y),function(h) S_pred(y[h],xx,tr))
  eval.pred.mean = unlist(lapply(eval.pred,mean))
  return(eval.pred.mean)
}
