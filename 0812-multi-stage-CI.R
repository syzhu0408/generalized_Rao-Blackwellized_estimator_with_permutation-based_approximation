##### drop-the-loser: normal K treatment J stage 1 control #####
# c(6,6,6)-60-2.296974
rm(list = ls())
# mu_treatment<- c(0,0,2)
# mu_control <- 0
# sigma_treatment <- c(6,6,6)
# sigma_control <- 6
# t <- c(1,2,3)/3
# unit_n <- 60
# J_stage <- c(3,2,1)
library(BSDA)
a <- function(sim,mu_treatment,sigma_treatment,sigma_control,t,unit_n,J_stage,alpha){
  set.seed(20240922+sim)
  J <- length(J_stage)
  c <- 2.296974
  mu_control <- 0
  K <- length(mu_treatment)
  n <- ceiling(t*length(t)*unit_n)
  all_da <- mvrnorm(n[J],mu=c(mu_control,mu_treatment),
                    Sigma = diag(c(sigma_control,sigma_treatment)**2,nrow = K+1,ncol = K+1))
  z_f <- function(da_treat,da_control,k) {z.test(da_treat,da_control,
                                                 sigma.x=sigma_treatment[k],sigma.y=sigma_control)$statistic}
  z <- mapply(function(k) mapply(function(j)  z_f(all_da[c(1:n[j]),k+1],all_da[c(1:n[j]),1],k),j=1:J) ,k=1:K)
  con <- matrix(FALSE,nrow = J, ncol = K)
  con[1,] <- TRUE
  for(j in 1:(J-1)){
    con[j+1,which(con[j,])[rank(z[j,con[j,]])>=(J_stage[j]-J_stage[j+1]+1)]] <- TRUE
  }
  end_stage <- apply(con, 2, sum)
  select_index <- which(end_stage==J)
  naive <- mean(all_da[,select_index+1])-mean(all_da[,1])
  stage_end <- mean(all_da[(n[J-1]+1):n[J],select_index+1])-mean(all_da[(n[J-1]+1):n[J],1])
  CI_f <- function(est,n,select_index) {        
    SE <- sqrt((sigma_treatment[select_index]**2+sigma_control**2)/n)
    c(est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
  }
  UMCVUE_per <- function(per_num){
    per_sample_treat <- list()
    for (k in 1:K) {
      per_sample_treat[[k]] <- replicate(per_num,sample(all_da[1:n[end_stage[k]],k+1],
                                                        n[end_stage[k]],FALSE), simplify = T)
    }
    per_sample_control <- replicate(per_num,sample(all_da[1:n[J],1],
                                                   n[J],FALSE), simplify = T)
    z_treat_I <- mapply(function(k) mapply(function(i) z_f(per_sample_treat[[k]][c(1:n[1]),i],
                                                           per_sample_control[c(1:n[1]),i],k),i=1:per_num),k=1:K)
    sim_index <- which(mapply(function(i) min(z_treat_I[i,con[2,]]) > max(z_treat_I[i,!con[2,]]) , i=1:per_num))
    re <- rep(NA,4)
    if(length(sim_index) > 1){
      z_treat_II <- mapply(function(k) mapply(function(i) z_f(per_sample_treat[[k]][c(1:n[2]),i],
                                                              per_sample_control[c(1:n[2]),i],k),i=sim_index),
                           k = which(con[2,]))
      sim_index2 <-  sim_index[mapply(function(i) min(z_treat_II[i,con[3,][con[2,]]]) > max(z_treat_II[i,!(con[3,][con[2,]])]) ,
                                      i=1:length(sim_index))]
      if(length(sim_index2) > 1){
        accept <- apply(per_sample_treat[[select_index]][c((n[J-1]+1):n[J]),sim_index2],2,mean)-
          apply(per_sample_control[c((n[J-1]+1):n[J]),sim_index2],2,mean)
        est <- mean(accept)
        SE <- sqrt((sigma_treatment[select_index]**2+sigma_control**2)/unit_n-var(accept))
        re <-  c(length(sim_index2)/per_num,est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
      }}
    re
  }
  UMCVUE_pe <- UMCVUE_per(per_num=10000)
  c(z[J,select_index]>c,end_stage,select_index,
    naive,CI_f(naive,n[J],select_index),stage_end,CI_f(stage_end,n[1],select_index),UMCVUE_pe)
}
set.seed(20240919)
mu_treat_matrix <- matrix(c(c(0.8,1.5,2.6),c(0,0,2)),3)
for(s in 1:dim(mu_treat_matrix)[2]){
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
                    ac(i,mu_treatment=mu_treat_matrix[,s],sigma_treatment=rep(6,3),sigma_control=6,
                       t=c(1,2,3)/3,unit_n=60,J_stage=c(3,2,1),alpha=0.05)
                  }
  close(pb)
  stopCluster(cls)
  time2 <- Sys.time()
  print(time2-time1)
  setwd("E:/0621/桌面/0919/new-res/multi stage")
  write.csv(t(res),paste("normal-","setting-",s,".csv",sep = ""))
}

##### MAMS group sequential desgin: binomal K treatment J stage 1 control #####
rm(list = ls())
# library(MAMS)
# design <- ordinal.mams(prob=c(0.02,0.98), or=12.25, or0=1, K=3, J=3, alpha=0.025,
#                  power=0.8, r=1:3, r0=1:3, ushape="obf",
#                  lshape="obf")
# low_bround <- design$l
# up_bround <- design$u
# unit_n <- design$n
a <- function(sim,p_treat,p_control,t,unit_n,low_bround,up_bround,alpha){
  set.seed(20240922+sim)
  K <- length(p_treat)
  n <- ceiling(t*length(t)*unit_n)
  J <- length(t)
  z_f <- function(da_treat,da_control,n_treat,n_control) {
    da_treat[which(da_treat==0)] <- 0.5
    da_control[which(da_control==0)] <- 0.5
    da_treat[which(da_treat==n_treat)] <- n_treat-0.5
    da_control[which(da_control==n_control)] <- n_control-0.5
    p_pool <- (da_treat+da_control)/(n_treat+n_control)
    (da_treat/n_treat-da_control/n_control)/
      sqrt(((p_pool)*(1-p_pool)*(1/n_treat+1/n_control)))
  }
  all_da <- apply(mapply(function(x) mapply(function(i) rbinom(1,c(n[1],diff(n))[x],
                                                               c(p_control,p_treat)[i]),
                                            i=1:(K+1)), x=1:J),1,cumsum)
  z <- t(mapply(function(j) z_f(da_treat=all_da[j,-1],da_control=all_da[j,1],
                                n_treat=n[j],n_control=n[j]), j=1:J))
  bround_f <- function(x,low,up){
    re <- rep(0,length(x))
    re[which(x < low)] <- -1
    re[which(x > up)] <- 1
    re
  }
  decide_matrix <- t(mapply(function(j) bround_f(z[j,],low=low_bround[j],up=up_bround[j]), j=1:J))
  final <- apply(decide_matrix,2,function(x) c(which(x!=0)[1],x[which(x!=0)[1]]))
  if(any(final[2,]==1)){
    end_stage <- min(final[1,final[2,]==1])
    recom_index <- which.max(z[end_stage,])
  }else if(all(final[2,]==-1)){
    end_stage <- max(final[1,])
    recom_index <- 0
  }
  final[2,final[1,] > end_stage] <- 0
  final[1,final[1,] > end_stage] <- end_stage
  
  naive <- mapply(function(k) all_da[final[1,k],k+1]/n[final[1,k]]-all_da[end_stage,1]/n[end_stage], k=1:K)
  stage1 <-  (all_da[1,-1]-all_da[1,1])/n[1]
  p_bar <- mapply(function(k) mean(c(all_da[final[1,k],k+1]/n[final[1,k]],all_da[end_stage,1]/n[end_stage])), k=1:K)
  CI_f <- function(est,N,index){
    SE <- sqrt(2*p_bar[index]*(1-p_bar[index])/N)
    c(est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
  }
  UMCVUE_per <- function(per_num){
    per_sample_treat <- list()
    for(k in 1:K){
      stage <- final[1,k]
      per_t <- replicate(per_num,sample(c(rep(1,all_da[stage,k+1]),
                                          rep(0,n[stage]-all_da[stage,k+1])),
                                        n[stage],FALSE), simplify = T)
      per_sample_treat[[k]] <- mapply(function(j) apply(per_t[1:n[j],], 2, sum),
                                      j=1:stage)
    }
    per_c <- replicate(per_num,sample(c(rep(1,all_da[end_stage,1]),
                                        rep(0,n[end_stage]-all_da[end_stage,1])),
                                      n[end_stage],FALSE), simplify = T)
    per_sample_control <- mapply(function(j) apply(per_c[1:n[j],], 2, sum),
                                 j=1:end_stage)
    
    index_f <- function(sample_treat,sample_control,stage,k){
      z_per <- t(mapply(function(j) z_f(da_treat=sample_treat[,j],da_control=sample_control[,j],
                                        n_treat=n[j],n_control=n[j]), j=1:stage))
      apply(z_per,2,function(x) identical(between(x,low_bround[1:stage],
                                                  up_bround[1:stage]),
                                          between(z[1:stage,k],low_bround[1:stage],
                                                  up_bround[1:stage])))
    }
    # sample_treat <- per_sample_treat[[k]]; sample_control <- per_sample_control;stage <- final[1,k]
    sim_index <- apply(sapply(X=1:K, function(k) index_f(per_sample_treat[[k]],
                                                         per_sample_control,final[1,k],k)),
                       1, all)
    mean_per <- sapply(X=1:K, function(k) (per_sample_treat[[k]][,1]-per_sample_control[,1])/n[1])
    
    if(!is.na(mean(sim_index))){
      accpet <- mean_per[sim_index,]
      est <- apply(mean_per[sim_index,],2,mean)
      SE <- sqrt(2*p_bar*(1-p_bar)/n[1]-apply(mean_per[sim_index,],2,var))
      c(mean(sim_index),est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
    }else{ rep(NA,10) }
  }
  if(end_stage==1){
    UMCVUE_pe <- c(NA,naive,as.numeric(mapply(index=1:K,function(index) CI_f(naive[index],n[final[1,index]],index)))[c(1,3,5)],
                   as.numeric(mapply(index=1:K,function(index) CI_f(naive[index],n[final[1,index]],index)))[c(2,4,6)])
  }else{
    UMCVUE_pe <- UMCVUE_per(per_num=10000)
  }
  c(end_stage,recom_index,final[2,],c(naive,as.numeric(mapply(index=1:K,function(index) CI_f(naive[index],n[final[1,index]],index))),
                                      stage1,as.numeric(mapply(index=1:K,function(index) CI_f(stage1[index],n[1],index))),
                                      UMCVUE_pe))
}
ac <- cmpfun(a)
set.seed(20240919)
p_treat_matrix <- matrix(c(c(0.02,0.10,0.20),c(0.02,0.02,0.20)),3)
for(s in 1:dim(p_treat_matrix)[2]){
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
                    ac(i,p_treat=p_treat_matrix[,s],p_control=0.02,
                       t=c(1,2,3)/3,unit_n=56,
                       low_bround=c(-4.115443, -2.910058, 2.376052),
                       up_bround=c(4.115443, 2.910058, 2.376052),alpha=0.05)
                  }
  close(pb)
  stopCluster(cls)
  time2 <- Sys.time()
  print(time2-time1)
  setwd("E:/0621/桌面/0919/new-res/multi stage")
  write.csv(t(res),paste("binary-","setting-",s,".csv",sep = ""))
  
}

##### select-the-best: gamma K treatment J stage 1 control #####
## sit3--stage
# c(2,2,2,2)--25
# c(2,3,3,3)--6
# (2,4,4,4)--4
# (2,1,1.5,2)--33
# (1,1,1,1)--55
# (2,1,2,3)--15
rm(list = ls())
t <- c(1,2,3,4)/4
J <- length(t)
library(rpact)
low_bround <- up_bround <- numeric(length = J)
d <- getDesignGroupSequential(kMax = length(t), informationRates = t, 
                              typeOfDesign = "asOF",typeBetaSpending = "bsOF", 
                              bindingFutility = T,alpha = 0.025, beta = 0.2, 
                              sided = 1)
up_bround <- d$criticalValues
low_bround <- c(d$futilityBounds,d$criticalValues[J])
up_bround[1] <- Inf
low_bround[1] <- -Inf
a <- function(sim,shape_treat,shape_control,scale_treat,scale_control,unit_n,low_bround,up_bround,alpha){
  set.seed(20240922+sim)
  library(BSDA)
  K <- length(scale_treat)
  n <- ceiling(t*length(t)*unit_n)
  w_p <- function(da_treat, da_control) wilcox.test(da_treat, da_control, alternative = "greater")$p.value
  p_adj_num <- function(num,p_v){
    com_ind <- combn(1:K,num)
    in_f <- function(x) select_index %in% x
    p_adj <- function(ind,p_v) min(p.adjust(p_v[ind],method="hochberg"))
    ind <- com_ind[,apply(com_ind,2,in_f)]
    if(is.null(dim(ind))){p_adj(ind, p_v=p_v)
    }else{max(apply(ind, 2, p_adj, p_v=p_v))}
  }
  all_da <- mapply(function(i) rgamma(n[J], shape=c(shape_control,shape_treat)[i], 
                                      scale = c(scale_control,scale_treat)[i]),i=1:(K+1))
  select_index <- which.max(apply(all_da[c(1:n[1]),],2,mean)/c(shape_control,shape_treat))-1
  if(select_index==0){end_stage <- 1}else{
    p_stageI <- mapply(function(k) w_p(all_da[c(1:n[1]),k+1],all_da[c(1:n[1]),1]),k=1:K)
    p_select_adjI <- max(sapply(X=1:K,p_adj_num,p_v=p_stageI)) 
    p_select <- mapply(function(j) w_p(all_da[c((n[j-1]+1):n[j]),select_index+1],
                                       all_da[c((n[j-1]+1):n[j]),1]),j=2:J)
    w <- sqrt(c(n[1],diff(n)))
    z_select <- c(qnorm(1-p_select_adjI),
                  mapply(function(j) sum((w*c(qnorm(1-p_select_adjI),qnorm(1-p_select)))[1:j])/
                           sqrt(sum(w[1:j]^2)),j=2:J))
    end_stage <- first(which(between(z_select,low_bround,up_bround)==F))
    if(z_select[end_stage] >= up_bround[end_stage]) {decide="E"} else {decide="F"}
  }
  if(end_stage == 1){
    c(end_stage,select_index,rep(NA,4))
  }else{
    select_da <- all_da[c(1:n[end_stage]),select_index+1]
    control_da <- all_da[c(1:n[end_stage]),1]
    naive <- mean(select_da)/shape_treat[select_index]-mean(control_da)/shape_control
    stage2 <- mean(all_da[c((n[1]+1):n[2]),select_index+1])/shape_treat[select_index]-
      mean(all_da[c((n[1]+1):n[2]),1])/shape_control
    CI_f <- function(est,n,select_index) {        
      SE <- sqrt(mean(select_da)^2 / (n * shape_treat[select_index]^3) + mean(control_da)^2 / (n * shape_control^3))
      c(est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
    }
    UMCVUE_per2 <- function(sim_sample_treat,sim_sample_control,sim_num,end_stage){
      sim_m_stageI <-  matrix(apply(all_da[c(1:n[1]),],2,mean),nrow = K+1,ncol = sim_num,byrow = F)
      sim_m_stageI[1,] <- apply(sim_sample_control[1:n[1],],2,mean)
      sim_m_stageI[select_index+1,] <- apply(sim_sample_treat[1:n[1],],2,mean)
      sim_index <- apply(sim_m_stageI/c(shape_control,shape_treat), 2, which.max)==(select_index+1)
      sim_index2 <- NULL
      if(sum(sim_index) > 1){
        sim_p_stageI <-  mapply(function(k) mapply(function(i) w_p(all_da[c(1:n[1]),k+1],
                                                                   sim_sample_control[c(1:n[1]),i]),i=which(sim_index)), k=1:K)
        sim_p_stageI[,select_index] <- mapply(function(i) w_p(sim_sample_treat[c(1:n[1]),i],
                                                              sim_sample_control[c(1:n[1]),i]),i=which(sim_index))
        sim_p_select_adjI <- apply(sim_p_stageI,1,function(x) max(sapply(X=1:K,p_adj_num,p_v=x)))
        sim_p_select <- mapply(function(j)
          mapply(function(i) w_p(sim_sample_treat[c((n[j-1]+1):n[j]),i],
                                 sim_sample_control[c((n[j-1]+1):n[j]),i]),i=which(sim_index)),j=2:end_stage)
        sim_p_select <- cbind(sim_p_select_adjI,sim_p_select)
        sim_z_select <- rbind(qnorm(1-sim_p_select[,1]),apply(qnorm(1-sim_p_select),
                                                              1,function(x) mapply(function(j)
                                                                sum((w[1:j]*x[1:j]))/sqrt(sum(w[1:j]^2)),j=2:end_stage)))
        sim_index2 <- apply(sim_z_select[1:end_stage,],2,function(x) 
          identical(between(x,low_bround[1:end_stage],
                            up_bround[1:end_stage]),
                    c(rep(TRUE,end_stage-1),FALSE)))
      }
      if(sum(sim_index2) > 1){
        accept <- (apply(sim_sample_treat[c((n[1]+1):n[2]),sim_index],2,mean)/shape_treat[select_index]-
                     apply(sim_sample_control[c((n[1]+1):n[2]),sim_index],2,mean)/shape_control)[sim_index2]
        est <- mean(accept)
        SE <- sqrt(mean(select_da)^2 / (unit_n * shape_treat[select_index]^3) + mean(control_da)^2 / (unit_n * shape_control^3)-var(accept))
        c(sum(sim_index2)/sim_num,est,est+qnorm(alpha/2)*SE,est-qnorm(alpha/2)*SE)
      } else{rep(NA,4)}
    }
    per_num <- 10000
    per_sample_treat <- replicate(per_num,sample(select_da,n[end_stage],FALSE), simplify = T)
    per_sample_control <- replicate(per_num,sample(control_da,n[end_stage],FALSE), simplify = T)
    UMCVUE_pe <- UMCVUE_per2(per_sample_treat,per_sample_control,per_num,end_stage)
    # sim_sample_treat=per_sample_treat;sim_sample_control=per_sample_control;sim_num=per_num
    c(end_stage,select_index,c(naive,CI_f(naive,n[end_stage],select_index),
                               stage2,CI_f(stage2,unit_n,select_index),UMCVUE_pe))
  }
}
ac <- cmpfun(a)
set.seed(20240919)
scale_treat_matrix <- matrix(c(c(1,1.2,1.5),c(1,1,1.5)),3)
for(s in 1:dim(scale_treat_matrix)[2]){
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
                    ac(i,shape_treat=c(2,2,2),shape_control=2,
                       scale_treat=scale_treat_matrix[,s],scale_control=1,unit_n=25,
                       low_bround,up_bround,alpha=0.05)
                  }
  close(pb)
  stopCluster(cls)
  time2 <- Sys.time()
  print(time2-time1)
  setwd("E:/0621/桌面/0919/new-res/multi stage")
  write.csv(t(res),paste("gamma-","setting-",s,".csv",sep = ""))
}