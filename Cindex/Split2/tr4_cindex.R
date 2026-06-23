# Cindex - transition 4
# rm(list = ls())


# set directory work ------------------------------------------------------
setwd("..../Split2")


# functions ---------------------------------------------------------------
source("bsSizeLN_functions.R")


# data --------------------------------------------------------------------
load("Cindexinfo.Rdata")
tr=4

# elements for posterior mean
y.x = to.cindex.info[[tr]]$time.x
xx.x = cbind(rep(1,nrow(to.cindex.info[[tr]])),
             to.cindex.info[[tr]]$age.x,to.cindex.info[[tr]]$ln.x,
             to.cindex.info[[tr]]$bs1.x,to.cindex.info[[tr]]$bs2.x,to.cindex.info[[tr]]$bs3.x,
             to.cindex.info[[tr]]$gr2.x,to.cindex.info[[tr]]$gr3.x,
             to.cindex.info[[tr]]$er.x)

y.y = to.cindex.info[[tr]]$time.y
xx.y = cbind(rep(1,nrow(to.cindex.info[[tr]])),
             to.cindex.info[[tr]]$age.y,to.cindex.info[[tr]]$ln.y,
             to.cindex.info[[tr]]$bs1.y,to.cindex.info[[tr]]$bs2.y,to.cindex.info[[tr]]$bs3.y,
             to.cindex.info[[tr]]$gr2.y,to.cindex.info[[tr]]$gr3.y,
             to.cindex.info[[tr]]$er.y)

hp.x=to.cindex.info[[tr]]$hosp.x
hp.y=to.cindex.info[[tr]]$hosp.y


# computation for x and y -------------------------------------------------
# tr4_f.x=sapply(1:length(y.x),function(h)  {eval_p_mean(y.x[h],xx=xx.x[h,],tr,hp=hp.x[h])}) 
# tr4_f.y=sapply(1:length(y.y),function(h)  {eval_p_mean(y.y[h],xx=xx.y[h,],tr,hp=hp.y[h])}) 


# tr4_cindex = cbind(tr4_f.x,tr4_f.y)
tr4_cindex = sapply(1:nrow(to.cindex.info[[tr]]),
                    function(h) {  eval_p_cindex(y.x=y.x[h],xx.x=xx.x[h,],tr.x=tr,hp.x=hp.x[h],
                                                 y.y=y.y[h],xx.y=xx.y[h,],tr.y=tr,hp.y=hp.y[h]) })
# file --------------------------------------------------------------------
save(tr4_cindex,file=paste0("tr",tr,"_cindex.Rdata"))

