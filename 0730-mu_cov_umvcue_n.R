mu_cov_UMVCUE <- function(control_eff,treat_eff,control_sig=NULL,treat_sig=NULL,
                          select_index,end_stage,t,K,J,n,type){

  control <- control_sig; treat <- treat_sig
  if(type=="binary"){
    V <- function(eff,j,i) eff*(1-eff)/n[j,i]
  }else{
    V <- function(sig,j,i) sig**2/n[j,i]
  }
  I <- function(i,j){1/(V(treat[i],j,i+1)+V(control,j,1))} 
  mu <- rep(0,end_stage+2)
  cov <-cbind(matrix(nrow = end_stage+2,ncol = end_stage-1),
              matrix(sqrt(I(select_index,1)),nrow = end_stage+2,ncol = 1),
              matrix(sqrt(I(select_index,end_stage)),nrow = end_stage+2,ncol = 2))
  for(i in 1:(end_stage-1)){
    for(j in i:(end_stage-1)){
      cov[i,j] <- sqrt(I(select_index,i)/I(select_index,j))
    }
  }
  for(i in 1:(end_stage+2)){
    for(j in (end_stage:(end_stage+2))){
      if(i <= end_stage-1){cov[i,j] <- cov[i,j] * sqrt(I(select_index,i))}
      if(i == end_stage){cov[i,j] <- cov[i,j] * sqrt(I(select_index,1))}
      if(i >= end_stage+1){cov[i,j] <- cov[i,j] * sqrt(I(select_index,end_stage))}
    }
  }
  # cov <- matrix(1,4,4)
  for(i in 1:(end_stage-1)){
    for(j in (end_stage:(end_stage+2))){
      if(j == end_stage){cov[i,j] <- -cov[i,j] * V(control,i,1)}
      if(j == end_stage+1){cov[i,j] <- -cov[i,j] * V(control,end_stage,1)}
      if(j == end_stage+2){cov[i,j] <- cov[i,j] * V(treat[select_index],end_stage,select_index+1)}
    }
  }
  diag(cov[(end_stage):(end_stage+2),(end_stage):(end_stage+2)]) <- diag(
    cov[(end_stage):(end_stage+2),(end_stage):(end_stage+2)])*
    c(V(control,1,1),V(control,end_stage,1),V(treat[select_index],end_stage,select_index+1))
  cov[end_stage,(end_stage+1):(end_stage+2)] <- cov[end_stage,(end_stage+1):(end_stage+2)]*
    c(V(control,end_stage,1),0)
  cov[end_stage+1,end_stage+2] <- 0
  cov[lower.tri(cov)] <- t(cov)[lower.tri(cov)]
  list(mu=mu,cov=cov)
}
