rm(list = ls())

# libreries ---------------------------------------------------------------
library(survival)
library(ggplot2)
#library(survminer)
library(splines)

library(rjags)
library(R2jags)


# data --------------------------------------------------------------------


load("res_up_II.RData")
res <- res
head(res)
res$time <- res$Tstop-res$Tstart


#--------------------------------------------------------------------------
#In the article Dynamics of breast-cancer relapse reveal late-recurring ER-positive genomic 
#subgroups, additional transitions were considered by stratifying them according to ER status 
#in order to avoid violating the Cox proportional hazards assumption. 
#In our case, we retained the original 9 transitions defined by the diagram.
#---------------------------------------------------------------------------
res2 <- res
res2$trans2 <- res2$trans
res2$trans2[res2$trans2==10] <- 1;res2$trans2[res2$trans2==11] <- 2
res2$trans2[res2$trans2==12] <- 3;res2$trans2[res2$trans2==13] <- 4
res2$trans2[res2$trans2==14] <- 5;res2$trans2[res2$trans2==15] <- 6
res2$trans2[res2$trans2==16] <- 7;res2$trans2[res2$trans2==17] <- 8
res2$trans2[res2$trans2==18] <- 9

res2$ER.status <-ifelse(res2$ER=='ER+',1,0)

res2.na <- na.omit(res2[,c(1:11,14:17)])


library(tidyverse)



X <- model.matrix(~-1+AGE+LN_O+bs(Size,df=3)+as.factor(Grade)+ER.status, data = res2.na)

N <- length(unique(res2.na$id))
ind.test <- unique(res2.na$id)
Ni <- unlist(lapply(1:N, function(i,...){
  length(res2.na$id[res2.na$id==ind.test[i]])
}))

age <- matrix(NA,N,max(Ni))
ln <- matrix(NA,N,max(Ni))
er <- matrix(NA,N,max(Ni))
size <- matrix(NA, N ,max(Ni))
grade <- matrix(NA, N ,max(Ni))
bs.size1 <- matrix(NA,N,max(Ni))
bs.size2 <- matrix(NA,N,max(Ni))
bs.size3 <- matrix(NA,N,max(Ni))
time.test <- matrix(NA,N,max(Ni))
cen.test <- matrix(NA,N,max(Ni))
trans.test <- matrix(NA,N,max(Ni))
hosp.test <- matrix(NA,N,max(Ni))


## solo para crear is.censored
time.is.new<- res2.na$time
cens <- res2.na$time
time.is.new[res2.na$status==0] <- NA
is.censored.new <- as.numeric(is.na(time.is.new))
is.censored.test <- matrix(NA,N,max(Ni))

count <- 1
for(i in 1:N){
  age[i,1:Ni[i]] <- res2.na[count:(Ni[i]+count-1),"AGE"] 
  ln[i,1:Ni[i]] <- res2.na[count:(Ni[i]+count-1),"LN_O"] 
  er[i,1:Ni[i]] <- res2.na[count:(Ni[i]+count-1),"ER.status"]
  size[i,1:Ni[i]] <- res2.na[count:(Ni[i]+count-1),"Size"]
  grade[i,1:Ni[i]] <- res2.na[count:(Ni[i]+count-1),"Grade"]
  bs.size1[i,1:Ni[i]] <- X[count:(Ni[i]+count-1),3]
  bs.size2[i,1:Ni[i]] <- X[count:(Ni[i]+count-1),4]
  bs.size3[i,1:Ni[i]] <- X[count:(Ni[i]+count-1),5]
  time.test[i,1:Ni[i]] <- time.is.new[count:(Ni[i]+count-1)] 
  cen.test[i,1:Ni[i]] <- cens[count:(Ni[i]+count-1)]  
  trans.test[i,1:Ni[i]] <- res2.na[count:(Ni[i]+count-1),"trans2"] 
  hosp.test[i,1:Ni[i]] <- res2.na[count:(Ni[i]+count-1),"hospital"] 
  is.censored.test[i,1:Ni[i]] <- is.censored.new[count:(Ni[i]+count-1)] 
  count <- count + Ni[i]
}

transition = matrix(NA,ncol = 9,nrow = nrow(trans.test))
for(i in 1:9){
  transition[,i] = apply(trans.test,1,function(h) any(h==i,na.rm = T))  
}
transition=as.data.frame(transition)
names(transition) <- paste0("trans_",1:9)

XXL <- model.matrix(~age[,1]+ln[,1]+bs.size1[,1]+bs.size2[,1]+bs.size3[,1]+as.factor(grade[,1])+er[,1])


##bnp model ---------------------------------------------------------------

bnp_surv_model <- "
model {
## stick
  for(g in 1:G){
       B[1,g] <- b[1,g]
  }

 for(h in 2:H){
  for(g in 1:(G-1)){
    B[h,g] <- b[h,g]*((1- b[h-1,g])*B[h-1,g]/b[h-1,g])
  }
  B[h,G] <- 1 - sum(B[h,1:(G-1)]) 
 }

## prob for selection index
  for (h in 1:H){
    for(g in 1:G){
       b[h,g] ~ dbeta(1, M)          
       BetaProb[h,g] <- B[h,g]/(sum(B[,g]))
    }
  }

## atoms 
  for (h in 1:H){
    for(g in 1:G){
      betah[1:p,g,h] ~ dmnorm(mub,Sigmab) 
    }
  }

## group index selection
  for (l in 1:n){ 
    for(j in 1:ID[l]){
     m[l,j] ~ dcat(BetaProb[,transition[l,j]])
    }
  }  
      
## likelihood   
for (l in 1:n) {
  for(j in 1:ID[l]){
    is.censored[l,j] ~ dinterval(time[l,j],cen[l,j])
    time[l,j] ~ dlnorm(mu[l,j],tau2[transition[l,j]])
    mu[l,j] <- inprod(betah[,transition[l,j],m[l,j]],X[l,]) + bk[hospital[l]] 
  }
}  

## prior for random effect hospital variance
  for(k in 1:K){
    sigma2k[k] <- 1/tau2k[k]
    tau2k[k] ~ dgamma(tk1[k],tk2[k])
     bk[k] ~ dnorm(0, tau2k[k])
  }
 

#hyperparameters
M ~ dgamma(a0,b0)
mub ~ dmnorm(m0,S0)
Sigmab ~ dwish(Psiinv,psi)


## Priors for sigma according the transition
 for(g in 1:G){
  sigma2[g] <- 1/tau2[g]
  tau2[g] ~ dgamma(gamag[g],deltag[g])
 }
}
"


# especification ----------------------------------------------------------
H0=25
#zellner's prior
g0 <- 9*N
Sigma0 <- g0*solve(t(XXL)%*%XXL)
beta0 <- rep(0,9)

d.bnp2 <- list(n = N,time=time.test,cen =cen.test, ID = Ni, 
               X = XXL, is.censored = is.censored.test, 
               transition = trans.test, hospital = hosp.test[,1],
               p=ncol(XXL),H=H0,G=9, K=5,
               a0 = 20, b0= 10, 
               m0 = beta0,S0 = Sigma0,
               Psiinv = Sigma0,psi=ncol(XXL)+2,
               gamag=rep(10,9),#rep(100,200,length.out=9),
               deltag=rep(2,9),#rep(100,200,length.out=9),
               tk1 = rep(12,9),#seq(100,200,length.out=5),
               tk2 = rep(12,9)#seq(100,200,length.out=5)
               )


p.bnp2 <- c("betah","M","sigma2","BetaProb","m","bk", "sigma2k")


out.DDP <- jags(data=d.bnp2,
                 parameters=p.bnp2,
                 model = textConnection(bnp_surv_model),
                 n.chains=2,
                 n.iter = 200000,
                 n.thin=10,
                 n.burnin=100000)   

save(out.DDP, file = paste0("Results_BNP_",H0,"_zellner.RData"))

