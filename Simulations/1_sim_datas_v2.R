#simulation data: 
rm(list = ls())
# functions - load --------------------------------------------------------
source("0_sim_functions_v2.R")

# simulation structure ----------------------------------------------------
ns <- c(250,500,1000)
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
cens.par = 2

# data: n1-par1 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[1],betas[1,],rates[1,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

# data: n1-par2 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[1],betas[2,],rates[2,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))


# data: n1-par3 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[1],betas[3,],rates[3,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))


# data: n2-par1 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[2],betas[1,],rates[1,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

# data: n2-par2 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[2],betas[2,],rates[2,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

# data: n2-par3 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[2],betas[3,],rates[3,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))


# data: n3-par1 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[3],betas[1,],rates[1,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

# data: n3-par2 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[3],betas[2,],rates[2,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

# data: n3-par3 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[3],betas[3,],rates[3,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))



# data: n1-par4 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[1],betas[4,],rates[4,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

# data: n2-par4 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[2],betas[4,],rates[4,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

# data: n3-par4 -----------------------------------------------------------
data_trans = list()
for(i in 1:100){
  data_trans[[i]] <- sim_data(ns[3],betas[4,],rates[4,],cens.par)
}
datas <- lapply(1:100,function(h) eval_X(h,datas = data_trans[[h]]))

