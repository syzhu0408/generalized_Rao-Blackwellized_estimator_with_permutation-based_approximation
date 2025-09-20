rm(list = ls())
##### gamma select-the-best #####
setwd("E:/0621/桌面/0919/new-res/multi")
for(j in c(1:2)){
  da <- read.csv(paste("gamma-","setting-",j,".csv",sep = ""))[,-1]
  colnames(da) <- c("end_stage","select_index","naive","naive_L","naive_U",
                    "stage-wise","stage2_L","stage2_U",
                    "permutation-index","permutation","permutation-L","permutation-U")
  name <- c("naive","stage-wise","permutation")
  theta <- scale_tdaat_matrix <- matrix(c(c(1,1.2,1.5),c(1,1,1.5)),3)-1
  
  RMSE <- bias <- CP <- Width <- matrix(ncol = length(name), nrow = 3)
  colnames(RMSE) <- colnames(bias) <- colnames(CP) <- colnames(Width) <- name
  for(i in 1:3){
    bias[i,] <- apply(da[da$select_index==i,name], 2, mean, na.rm=T)-theta[i,j]
    RMSE[i,] <- apply(da[da$select_index==i,name], 2, function(x) sqrt(mean((x-theta[i,j])^2, na.rm=T)))
    CP[i,] <- c(mean(apply(cbind(da[da$select_index==i,"naive_L"] < theta[i,j] , da[da$select_index==i,"naive_U"] > theta[i,j]),1,all),na.rm=T),
                mean(apply(cbind(da[da$select_index==i,"stage2_L"] < theta[i,j] , da[da$select_index==i,"stage2_U"] > theta[i,j]),1,all),na.rm=T),
                mean(apply(cbind(da[da$select_index==i,"permutation-L"] < theta[i,j] , da[da$select_index==i,"permutation-U"] > theta[i,j]),1,all),na.rm=T))
    Width[i,] <- c(mean(da[da$select_index==i,"naive_U"] - da[da$select_index==i,"naive_L"] ,na.rm=T),
                   mean(da[da$select_index==i,"stage2_U"] - da[da$select_index==i,"stage2_L"] ,na.rm=T),
                   mean(da[da$select_index==i,"permutation-U"] - da[da$select_index==i,"permutation-L"] ,na.rm=T))
  }
  bias <- as.data.frame(bias)
  RMSE <- as.data.frame(RMSE)
  CP <- as.data.frame(CP)
  Width <- as.data.frame(Width)
  write.csv(bias,paste(paste("gamma-bias",j,sep="-"),".csv",sep=""))
  write.csv(RMSE,paste(paste("gamma-RMSE",j,sep="-"),".csv",sep=""))
  write.csv(CP,paste(paste("gamma-CP",j,sep="-"),".csv",sep=""))
  write.csv(Width,paste(paste("gamma-Width",j,sep="-"),".csv",sep=""))
}

## plot
setwd("E:/0621/桌面/0919/new-res/multi")
j <- 2
da_bias <- read.csv(paste(paste("gamma-bias",j,sep="-"),".csv",sep=""))
da_RMSE <- read.csv(paste(paste("gamma-RMSE",j,sep="-"),".csv",sep=""))
da_CP <- read.csv(paste(paste("gamma-CP",j,sep="-"),".csv",sep=""))
da_Width <- read.csv(paste(paste("gamma-Width",j,sep="-"),".csv",sep=""))
da_bias$measure <- "conditional bias"
da_RMSE$measure <- "conditional RMSE"
da_CP$measure <- "conditional CP"
da_Width$measure <- "conditional Width"
colnames(da_bias) <- colnames(da_RMSE) <- colnames(da_CP) <- colnames(da_Width) <- c("select_index","niave","stage-wise","permutation","measure")

da <- reshape2::melt(rbind(da_bias,da_RMSE,da_CP,da_Width),id.vars=c("select_index","measure"),
                     c("niave","stage-wise","permutation"),
                     value.name="value",variable="method")
da$measure <- factor(da$measure,levels = c("conditional bias","conditional RMSE","conditional CP","conditional Width"))
da$treatment <- factor(da$select_index)
ref <- data.frame(measure = "conditional CP", yintercept = 0.95)
ref$measure <- factor(ref$measure,levels = c("conditional bias","conditional RMSE","conditional CP","conditional Width"))
pdf(file=paste("gamma-",j,".pdf",sep=""),width=9,height=9)
ggplot(da, aes(treatment,value,colour = method,group=method)) + 
  geom_point(size = 2) +   
  geom_line(linewidth = 1,aes(linetype = method)) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "bottom",
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold")
  )+
  ylab(NULL)+
  geom_hline(data = ref, aes(yintercept = yintercept), linetype = "dashed", color = "black")+
  facet_wrap(measure~.,scales="free",nrow = 2)
dev.off()


##### binomal MAMS GSD #####
setwd("E:/0621/桌面/0919/new-res/multi")
K <- 3
for(j in c(1:2)){
  da <- read.csv(paste("binary-","setting-",j,".csv",sep = ""))[,-1]
  names(da) <- c("end_stage","select_index",paste("index-",1:K,sep = ""),
                    paste("naive",1:K,sep = "."),as.vector(sapply(X=1:K, function(k) paste(c("naive.L","naive.U"),k,sep="."))) ,
                    paste("stage.wise",1:K,sep = "."),as.vector(sapply(X=1:K, function(k) paste(c("stage.wise.L","stage.wise.U"),k,sep="."))),
                    "per.RB.index",paste("per.RB",1:K,sep = "."), paste("per.RB.L",1:K,sep="."),paste("per.RB.U",1:K,sep="."))
  name <- c("naive","stage.wise","per.RB")
  theta <- matrix(c(c(0.02,0.10,0.20),c(0.02,0.02,0.20)),3)-0.02

  RMSE <- bias <- CP <- Width <- matrix(ncol = length(name), nrow = 3)
  colnames(RMSE) <- colnames(bias) <- colnames(CP) <- colnames(Width) <- name

  bias <- t(mapply(function(i) apply(da[, paste(name,i,sep = ".")], 2, mean, na.rm=T)-theta[i,j], i=1:K))
  RMSE <- t(mapply(function(i) apply(da[, paste(name,i,sep = ".")], 2, function(x) sqrt(mean((x-theta[i,j])^2, na.rm=T))), i=1:K))
  CP <- t(mapply(function(i) mapply(function(m) mean(apply(cbind(da[,paste(name[m],"L",i,sep = ".")] < theta[i,j], 
                                                                 da[,paste(name[m],"U",i,sep = ".")] > theta[i,j]),1,all),na.rm=T), 
                                    m = 1:length(name)),i=1:K))
  Width <-  t(mapply(function(i) mapply(function(m) mean(da[,paste(name[m],"U",i,sep = ".")]-da[,paste(name[m],"L",i,sep = ".")],na.rm=T), 
                                        m = 1:length(name)),i=1:K))
  bias <- as.data.frame(bias)
  RMSE <- as.data.frame(RMSE)
  CP <- as.data.frame(CP)
  Width <- as.data.frame(Width)
  write.csv(bias,paste(paste("binary-bias",j,sep="-"),".csv",sep=""))
  write.csv(RMSE,paste(paste("binary-RMSE",j,sep="-"),".csv",sep=""))
  write.csv(CP,paste(paste("binary-CP",j,sep="-"),".csv",sep=""))
  write.csv(Width,paste(paste("binary-Width",j,sep="-"),".csv",sep=""))
}

## plot
setwd("E:/0621/桌面/0919/new-res/multi")
j <- 1
da_bias <- read.csv(paste(paste("binary-bias",j,sep="-"),".csv",sep=""))
da_RMSE <- read.csv(paste(paste("binary-RMSE",j,sep="-"),".csv",sep=""))
da_CP <- read.csv(paste(paste("binary-CP",j,sep="-"),".csv",sep=""))
da_Width <- read.csv(paste(paste("binary-Width",j,sep="-"),".csv",sep=""))
da_bias$measure <- "bias"
da_RMSE$measure <- "RMSE"
da_CP$measure <- "CP"
da_Width$measure <- "Width"
colnames(da_bias) <- colnames(da_RMSE) <- colnames(da_CP) <- colnames(da_Width) <- c("select_index","niave","stage-wise","permutation","measure")

da <- reshape2::melt(rbind(da_bias,da_RMSE,da_CP,da_Width),id.vars=c("select_index","measure"),
                     c("niave","stage-wise","permutation"),
                     value.name="value",variable="method")
da$measure <- factor(da$measure,levels = c("bias","RMSE","CP","Width"))
da$treatment <- factor(da$select_index)
ref <- data.frame(measure = "CP", yintercept = 0.95)
ref$measure <- factor(ref$measure,levels = c("bias","RMSE","CP","Width"))
pdf(file=paste("binary-",j,".pdf",sep=""),width=9,height=9)
ggplot(da, aes(treatment,value,colour = method,group=method)) + 
  geom_point(size = 2) +   
  geom_line(linewidth = 1,aes(linetype = method)) +
  theme_bw()+
  theme(
    axis.title.x = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.key.size = unit(0.8, "cm"),
    legend.position = "bottom",
    strip.text.x = element_text(size = 12, face = "bold"),
    strip.text.y = element_text(size = 12, face = "bold")
  )+
  ylab(NULL)+
  geom_hline(data = ref, aes(yintercept = yintercept), linetype = "dashed", color = "black")+
  facet_wrap(measure~.,scales="free",nrow = 2)
dev.off()

##### normal drop-the-loser #####
setwd("E:/0621/桌面/0919/new-res/multi")
K <- 3
for(j in c(1:2)){
  da <- read.csv(paste("normal-","setting-",j,".csv",sep = ""))[,-c(1:2)]
  per_matrix <- gtools::permutations(n = 3, r = 3, v = c(1, 2, 3))
  names(da) <- c(paste("end_stage",1:K,sep = "."),"select_index",
                 "naive","naive.L","naive.U","stage.wise","stage.wise.L","stage.wise.U" ,
                 "per.RB.index","per.RB","per.RB.L","per.RB.U")
  name <- c("naive","stage.wise","per.RB")
  theta <- matrix(c(c(0.8,1.5,2.6),c(0,0,2)),3)
  
  RMSE <- bias <- CP <- Width <- matrix(ncol = length(name), nrow = dim(per_matrix)[1])
  colnames(RMSE) <- colnames(bias) <- colnames(CP) <- colnames(Width) <- name
  prob <- as.numeric()
  for(i in 1:dim(per_matrix)[1]){
    index <- apply(da[,paste("end_stage",1:K,sep = ".")], 1, function(x) identical(as.numeric(x), per_matrix[i,]))
    bias[i,] <- apply(da[index, name], 2, mean, na.rm=T)-theta[which.max(per_matrix[i,]),j]
    RMSE[i,] <- apply(da[index, name], 2, function(x) sqrt(mean((x-theta[which.max(per_matrix[i,]),j])^2, na.rm=T)))
    CP[i,] <- mapply(function(m) mean(apply(cbind(da[index,paste(name[m],"L",sep=".")] < theta[which.max(per_matrix[i,]),j], 
                                                  da[index,paste(name[m],"U",sep=".")] > theta[which.max(per_matrix[i,]),j]),1,all),na.rm=T), 
                     m = 1:length(name))
    Width[i,] <- mapply(function(m) mean(da[index,paste(name[m],"U",sep=".")]-da[index,paste(name[m],"L",sep=".")],na.rm=T), 
                        m = 1:length(name))
    prob <- c(prob,mean(index))
  }
  write.csv(round(cbind(per_matrix,prob,bias,RMSE,CP,Width),3),paste("normal-",j,".csv",sep=""))
}



