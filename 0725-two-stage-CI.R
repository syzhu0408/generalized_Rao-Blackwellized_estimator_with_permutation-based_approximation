library(MASS)
library(data.table)
library(compiler)
library(condMVNorm)
library(pbapply)
library(BSDA)
##### normal K treatment 2 stage 1 control #####
rm(list = ls())
#### sit1--sample size
## sit1--sample size
# 0.9--78
# 0.8--59
# 0.7--47
# 0.6--38
# 0.5--31
# 0.4--24
# 0.3--18
rm(list = ls())
a <- function(sim,mu_treatment,unit_n,sigma_treatment,sigma_control,alpha){
  set.seed(20240919+sim)
  K <- length(mu_treatment)
  mu_control <- 0
  n1 <- n2 <- unit_n
  z_p <- function(da_treat,da_control,k) {z.test(da_treat,da_control,alternative="greater",
                                                 sigma.x=sigma_treatment[k],sigma.y=sigma_control)$p.value}
  da_stage1 <- mvrnorm(n1,mu=c(mu_control,mu_treatment),
                       Sigma = diag(c(sigma_control,sigma_treatment)**2,nrow = K+1,ncol = K+1))
  p_v_stageI <- mapply(function(i) z_p(da_treat=da_stage1[,i+1],da_control=da_stage1[,1],i),i=1:K)
  select_index <- which.min(p_v_stageI)
  if(mean(da_stage1[,select_index+1])-mean(da_stage1[,1]) <= 0){ select_index <-0 }
  if(select_index==0){
    end_stage <- 1; c(end_stage,select_index,rep(NA,14))
  }else{
    end_stage <- 2
    da_stage2 <- mvrnorm(n2,mu=c(mu_control,mu_treatment[select_index]),
                         Sigma = diag(c(sigma_control,sigma_treatment[select_index])**2,nrow = 2,ncol = 2))
    select_da <- c(da_stage1[,select_index+1],da_stage2[,2])
    control_da <- c(da_stage1[,1],da_stage2[,1])
    
    naive <- mean(select_da)-mean(control_da)
    stage2 <- mean(da_stage2[,2])-mean(da_stage2[,1])
    CI_f <- function(est,n,select_index) {        
      SE <- sqrt((sigma_treatment[select_index]**2+sigma_control**2)/n)
      c(est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
    }
    
    UMVCUE_Stallard <- function(mc_num,select_index,control_sig,treat_sig){
      stage1 <- apply(da_stage1, 2, mean)
      source("E:/0621/桌面/12.18/0730-mu_cov_umvcue_n.R")
      I <- function(stage,index)  {
        if(stage==1){1/(control_sig**2/n1+treat_sig[index]**2/n1)}else{
          1/(control_sig**2/(n1+n2)+treat_sig[index]**2/(n1+n2))}
      }  
      r_select <- mean(select_da)*sqrt(I(2,select_index))
      r_control <- mean(control_da)*sqrt(I(2,select_index))
      
      mu_cov <- mu_cov_UMVCUE(control_eff=0,treat_eff=rep(0,K),
                              control_sig,treat_sig,
                              select_index,end_stage,
                              t=NA,K,J=2,
                              n=matrix(c(n1,n1+n2),nrow=2,ncol=K+1),
                              type="con")
      mc <- rcmvnorm(mc_num,mean=mu_cov$mu,sigma = mu_cov$cov,
                     dependent=1:end_stage,given=c((end_stage+1):(end_stage+2)),
                     X.given=c(r_control,r_select))
      mc_stage2 <- (sqrt(I(2,select_index))*(r_select-r_control)-sqrt(I(1,select_index))*mc[,1])/
        (I(2,select_index)-I(1,select_index))
      z_other <- mapply(function(i)(stage1[i+1]-mc[,2]/sqrt(I(1,select_index)))*sqrt(I(1,i)),
                        i=c(1:K)[-select_index])
      sim_index <- mc[,1] > apply(z_other, 1, max) & mc[,1] > 0
      if(sum(sim_index) > 1){
        accept <- mc_stage2[mc[,1] > apply(z_other, 1, max) & mc[,1] > 0]
        est <- mean(accept)
        SE <- sqrt((sigma_treatment[select_index]**2+sigma_control**2)/n2-var(accept))
        c(mean(sim_index),est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
      }else{
        rep(NA,4)
      }
    }
    sim_f <- function(sim_sample_treat,sim_sample_control,sim_num){
      sim_p_stageI <-  mapply(function(k) mapply(function(i) z_p(da_stage1[,k+1],
                                                                 sim_sample_control[c(1:n1),i],k),i=1:sim_num),
                              k=1:K)
      sim_p_stageI[,select_index] <- mapply(function(i) z_p(sim_sample_treat[c(1:n1),i],
                                                            sim_sample_control[c(1:n1),i],select_index),i=1:sim_num)
      sim_index <- (apply(sim_p_stageI, 1, which.min)==select_index) &
        (apply(sim_sample_treat[c(1:n1),],2,mean)-apply(sim_sample_control[c(1:n1),],2,mean) > 0)
      if(sum(sim_index) > 1){
        accept <- apply(sim_sample_treat[(n1+1):(n1+n2),sim_index], 2, mean)-apply(sim_sample_control[(n1+1):(n1+n2),sim_index], 2, mean)
        est <- mean(accept)
        SE <- sqrt((sigma_treatment[select_index]**2+sigma_control**2)/n2-var(accept))
        c(mean(sim_index),est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
        }else{
        rep(NA,4)
      }
    }
    sim_fc <- cmpfun(sim_f)
    per_num <- 10000
    UMCVUE_T <- UMVCUE_Stallard(mc_num=10000,select_index,control_sig = sigma_control, treat_sig = sigma_treatment)
    per_sample_treat <- replicate(per_num,sample(select_da,length(select_da),FALSE), simplify = T)
    per_sample_control <- replicate(per_num,sample(control_da,length(control_da),FALSE), simplify = T)
    # sim_sample_treat <- per_sample_treat; sim_sample_control <- per_sample_control;sim_num <- per_num
    UMCVUE_p <-sim_f(per_sample_treat,per_sample_control,per_num)
    c(end_stage,select_index,c(naive,CI_f(naive,n1+n2,select_index),stage2,CI_f(stage2,n2,select_index),UMCVUE_T,UMCVUE_p))
  }
}
set.seed(20240919)
for(unit_n in c(18,24,31,38,47,59,78)){
  sim_n <- 10
  
  time1 <- Sys.time()
  cls <- makeSOCKcluster(18)
  registerDoSNOW(cls)
  pb <- txtProgressBar(max=sim_n, style=3, char = "*",)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  ac <- cmpfun(a)
  res  <- foreach(i=c(1:sim_n),
                  .options.snow=opts,
                  .combine = cbind,
                  .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","BSDA")) %dopar% {
                    ac(i,mu_treatment=c(0.8,1.5,2.6),sigma_treatment=c(6,6,6),sigma_control=6,unit_n,alpha=0.05)
                  }
  close(pb)
  stopCluster(cls)
  time2 <- Sys.time()
  print(time2-time1)
  setwd("E:/0621/桌面/0919/new-res/2 stage")
  write.csv(t(res),paste("normal-",unit_n,".csv",sep = ""))
}

#### sit2--sigma
## sit2--sigma
# c(9,9,9)--95
# c(12,12,12)--145
# c(3,6,9)--88
# c(3,9,12)--131
# c(6,9,12)--134
rm(list = ls())
mu_treatment<- c(0.8,1.5,2.6)
mu_control <- 0
K <- length(mu_treatment)
sigma_treatment <- c(9,9,9)
sigma_control <- 6
unit_n <- 96
res <- pbsapply(X=c(1:5),ac)
write.csv(t(res),"C://Users//13379//Desktop//0919//sim//normal-var//sit2//c(9,9,9).csv")


##### binomal  K treatment 2 stage 1 control #####
## sit1--sample size
# 0.9--40
# 0.8--33
# 0.7--29
# 0.6--26
# 0.5--23
# 0.4--21
# 0.3--18
rm(list = ls())
a <- function(sim,p_treat,p_control,unit_n,alpha){
  set.seed(20240922+sim)
  n1 <- n2 <- unit_n
  K <- length(p_treat)
  stage1 <- mapply(function(i) rbinom(1,n1, c(p_control,p_treat)[i]),i=1:(K+1))
  sel_f <- function(x) min(which(x >= (max(x)-0.02)))-1
  select_index <- sel_f(stage1/n1)
  if(select_index==0){
    end_stage <- 1
    c(end_stage,select_index,rep(NA,14))
  }else{
    end_stage <- 2
    stage2 <- mapply(function(i) rbinom(1,n2, c(p_control,p_treat[select_index])[i]),i=1:2)
    all_da <- stage1[c(1,(select_index+1))]+stage2
    naive <- (all_da[2]-all_da[1])/(n1+n2)
    stage2 <- (stage2[2]-stage2[1])/n2
    p_bar <- mean(c(all_da[2]/(n1+n2),all_da[1]/(n1+n2)))
    CI_f <- function(est,N){
      SE <- sqrt(2*p_bar*(1-p_bar)/N)
      c(est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
    }
    mc_e <- function(m){
      mc_stage1 <- matrix(stage1,nrow = K+1, ncol = length(0:all_da[1])*length(0:all_da[2]), byrow = F)
      mc_stage1[1,] <- expand.grid(0:all_da[1],0:all_da[2])[,1]
      mc_stage1[select_index+1,] <- expand.grid(0:all_da[1],0:all_da[2])[,2]
      mc_index <- apply(mc_stage1/n1, 2, sel_f)==select_index
      prob <- dhyper(mc_stage1[select_index+1,mc_index],all_da[2],n1+n2-all_da[2],n1)*
        dhyper(mc_stage1[1,mc_index],all_da[1],n1+n2-all_da[1],n1)
      est <- sum(prob*((all_da[2]-mc_stage1[select_index+1,mc_index])-(all_da[1]-mc_stage1[1,mc_index])))/sum(prob)/n2
      SE <- sqrt(2*p_bar*(1-p_bar)/n2-(sum(prob/sum(prob)*((all_da[2]-mc_stage1[select_index+1,mc_index])/n2-
                                                                (all_da[1]-mc_stage1[1,mc_index])/n2)^2)-est^2))
      c(sum(prob),est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
    }
    per_f <- function(per_num){
      per_stage1 <- matrix(stage1,nrow = K+1, ncol = per_num, byrow = F)
      per_stage1[1,] <- rhyper(per_num,all_da[1],n1+n2-all_da[1],n1)
      per_stage1[select_index+1,] <- rhyper(per_num,all_da[2],n1+n2-all_da[2],n1)
      per_index <- apply(per_stage1/n1, 2, sel_f)==select_index
      if(sum(per_index) > 1){
        accept <- (all_da[2]-per_stage1[select_index+1,per_index])/n2-(all_da[1]-per_stage1[1,per_index])/n2
        est <- mean(accept)
        SE <- sqrt(2*p_bar*(1-p_bar)/n2-var(accept))
        c(mean(per_index),est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
      }else{
        rep(NA,4)
      }
    }
    c(end_stage,select_index,c(naive,CI_f(naive,n1+n2),stage2,CI_f(stage2,n2),
                               mc_e(mc_num),per_f(per_num=10000)))
  }
}
ac <- cmpfun(a)
sim_n <- 100000
set.seed(20240919)
for(unit_n in c(18,21,23,26,29,33,40)){
  sim_n <- 10
  
  time1 <- Sys.time()
  cls <- makeSOCKcluster(18)
  registerDoSNOW(cls)
  pb <- txtProgressBar(max=sim_n, style=3, char = "*",)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  ac <- cmpfun(a)
  res  <- foreach(i=c(1:sim_n),
                  .options.snow=opts,
                  .combine = cbind,
                  .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","BSDA")) %dopar% {
                    ac(i,p_treat=c(0.02,0.10,0.20),p_control=0.02,unit_n,alpha=0.05)
                  }
  close(pb)
  stopCluster(cls)
  time2 <- Sys.time()
  print(time2-time1)
  setwd("E:/0621/桌面/0919/new-res/2 stage")
  write.csv(t(res),paste("binary-",unit_n,".csv",sep = ""))
}

#### sit2--p_treatment
## sit2--p_treatment
# c(0.02,0.02,0.10)--91
# c(0.02,0.02,0.20)--32
# c(0.02,0.20,0.20)--26
# c(0.10,0.20,0.20)--27
# c(0.20,0.20,0.20)--23
rm(list = ls())
p_treat <- c(0.02,0.10,0.20)
p_control <- 0.02
n1 <- 18
n2 <- 18
sim_n <- 100000
per_num <- 10000
mc_num <- 10000
re <- pbsapply(X=c(1:sim_n),ac)
write.csv(t(re),"C://Users//13379//Desktop//0919//original//power//binomal-38.csv")

##### gamma K treatment 2 stage 1 control #####
#### sit1--samplesize
# 0.3--14
# 0.4--18
# 0.5--23
# 0.6--28
# 0.7--33
# 0.8--41
# 0.9--53
rm(list = ls())
a <- function(sim,shape_treat,shape_control,scale_treat,scale_control,unit_n,alpha){
  set.seed(20240922+sim)
  n1 <- n2 <- unit_n
  K <- length(shape_treat)
  stage1 <- mapply(function(i) rgamma(n1, shape=c(shape_control,shape_treat)[i], 
                                      scale = c(scale_control,scale_treat)[i]),i=1:(K+1))
  select_index <- which.max(apply(stage1,2,mean)/c(shape_control,shape_treat))-1
  if(select_index==0){
    end_stage <- 1
    c(end_stage,select_index,rep(NA,14))
  }else{
    end_stage <- 2
    stage2 <- mapply(function(i) rgamma(n2, shape=c(shape_control,shape_treat[select_index])[i], 
                                        scale = c(scale_control,scale_treat[select_index])[i]),i=1:2)
    select_da <- c(stage1[,select_index+1],stage2[,2])
    control_da <- c(stage1[,1],stage2[,1])
    naive_e <- mean(select_da)/shape_treat[select_index]-mean(control_da)/shape_control
    stage2_e <- mean(stage2[,2])/shape_treat[select_index]-mean(stage2[,1])/shape_control
    CI_f <- function(est,n,select_index) {        
      SE <- sqrt(mean(select_da)^2 / (n * shape_treat[select_index]^3) + mean(control_da)^2 / (n * shape_control^3))
      c(est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
    }
    mc <- function(N){
      mc_sample <- function(z,shape){
        f_f <- function(y) {((z-n2*y)/n1)^(n1*shape-1)*y^(n2*shape-1)}
        xy <- data.frame(proposed = runif(N, min = 0, max = ceiling(z/n2) ))
        xy$fit <- f_f(xy$proposed)
        xy$random <- runif(N, min = 0, max = 1)
        xy$accepted <- with(xy, random <= fit/max(xy$fit))
        xy[xy$accepted, ]$proposed
      }
      sel_stage2 <- mc_sample(sum(select_da),shape_treat[select_index])
      con_stage2 <- mc_sample(sum(control_da),shape_control)
      mc_n <- min(length(sel_stage2),length(con_stage2))
      index <- ((sum(select_da)-n2*sel_stage2[1:mc_n])/shape_treat[select_index] > 
                  (sum(control_da)-n2*con_stage2[1:mc_n])/shape_control) &
        ((sum(select_da)-n2*sel_stage2[1:mc_n])/n1/shape_treat[select_index] >
           max(apply(stage1[,-c(1,select_index+1)],2,mean)/shape_treat[-select_index]))
      accept <- sel_stage2[1:mc_n][index]/shape_treat[select_index]-con_stage2[1:mc_n][index]/shape_control
      est <- mean(accept)
      SE <- sqrt(mean(select_da)^2 / (n2 * shape_treat[select_index]^3) + mean(control_da)^2 / (n2 * shape_control^3)-var(accept))
      c(mean(index),est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
    }
    sim_f <- function(sim_sample_treat,sim_sample_control,sim_num){
      stage1_pre <- matrix(apply(stage1,2,mean),nrow = K+1,ncol = per_num)
      stage1_pre[1,] <- apply(sim_sample_control[1:n1,], 2, mean)
      stage1_pre[select_index+1,] <- apply(sim_sample_treat[1:n1,], 2, mean)
      per_index <- apply(stage1_pre/c(shape_control,shape_treat), 2, which.max)==(select_index+1)
      if(sum(per_index) > 1){
        accept <- apply(sim_sample_treat[(n1+1):(n1+n2),per_index], 2, mean)/shape_treat[select_index]-
          apply(sim_sample_control[(n1+1):(n1+n2),per_index], 2, mean)/shape_control
        est <- mean(accept)
        SE <- sqrt(mean(select_da)^2 / (n2 * shape_treat[select_index]^3) + mean(control_da)^2 / (n2 * shape_control^3)-var(accept))
        c(mean(per_index),est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
      }else{
        rep(NA,4)
      }
    }
    sim_fc <- cmpfun(sim_f)
    per_num <- 10000
    per_sample_treat <- replicate(per_num,sample(c(select_da),length(select_da),FALSE), simplify = T)
    per_sample_control <- replicate(per_num,sample(control_da,length(control_da),FALSE), simplify = T)
    # sim_sample_treat <- per_sample_treat; sim_sample_control <- per_sample_control;sim_num <- per_num
    UMCVUE_p <-sim_f(per_sample_treat,per_sample_control,per_num)
    mc_e <-mc(N = 1000000)
    c(end_stage,select_index,naive_e,CI_f(naive_e,n1+n2,select_index),
      stage2_e,CI_f(stage2_e,n1,select_index),mc_e,UMCVUE_p)
  }
}
ac <- cmpfun(a)
for(unit_n in c(14,18,23,28,33,41,53)){
  sim_n <- 10
  
  time1 <- Sys.time()
  cls <- makeSOCKcluster(18)
  registerDoSNOW(cls)
  pb <- txtProgressBar(max=sim_n, style=3, char = "*",)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  
  ac <- cmpfun(a)
  res  <- foreach(i=c(1:sim_n),
                  .options.snow=opts,
                  .combine = cbind,
                  .packages = c("MASS","data.table","compiler","condMVNorm","pbapply","BSDA")) %dopar% {
                    ac(i,shape_treat=c(2,2,2),shape_control=2,scale_treat=c(1,1.2,1.5),scale_control=1,unit_n,alpha=0.05)
                  }
  close(pb)
  stopCluster(cls)
  time2 <- Sys.time()
  print(time2-time1)
  setwd("E:/0621/桌面/0919/new-res/2 stage")
  write.csv(t(res),paste("gamma-",unit_n,".csv",sep = ""))
}

#### sit2--shape
## sit2--shape
# c(2,3,3,3)--9
# c(2,4,4,4)--5
# (2,1,1.5,2)--46
# (1,1,1,1)--90
# (2,1,2,3)--16
rm(list = ls())
shape_treat <- c(3,3,3)
shape_control <- 2
s_treat <- c(1,1.2,1.5)
s_control <- 1
unit_n <- 9
sim_n <- 100000
N <- 100000
per_num <- 10000
re <- pbsapply(X=c(1:sim_n),ac)
write.csv(t(re),"C:/Users/13379/Desktop/0919/sim/gamma//sit2//2-3-3-3.csv")


