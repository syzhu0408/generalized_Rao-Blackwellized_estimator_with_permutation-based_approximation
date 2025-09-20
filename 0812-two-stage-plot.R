rm(list = ls())
# label <- "normal"
label <- "binary"
# label <- "gamma"
all_bias <- all_RMSE <- all_CP <- all_Width <- data.frame()
setwd("E:/0621/桌面/0919/new-res/2 stage/")
for(i in 1:7){
  if(label=="normal"){unit_n <- c(18,24,31,38,47,59,78);theta <- c(0.8,1.5,2.6)
  }else if(label=="binary"){unit_n <- c(18,21,23,26,29,33,40);theta <- c(0.02,0.10,0.20)-0.02
  }else if(label=="gamma"){unit_n <- c(14,18,23,28,33,41,53);theta <- c(1,1.2,1.5)-1}
  re <- read.csv(paste(label,"-",unit_n[i],".csv",sep = ""))[,-1]
  
  colnames(re) <- c("end_stage","select_index",
                    "naive","naive_L","naive_U",
                    "stage2","stage2_L","stage2_U",
                    "RB-index","RB","RB-L","RB-U",
                    "per-RB-index","per-RB","per-RB-L","per-RB-U")
  name_CI <- name_point <- c("stage2","naive","RB","per-RB")
  CP <- Width <- RMSE <- bias <- matrix(ncol = length(name_point), nrow = 3)
  colnames(CP) <- colnames(Width) <- colnames(RMSE) <- colnames(bias) <- name_point
  
  for(i in 1:3){
    bias[i,] <- apply(re[re$select_index==i,name_point], 2, mean, na.rm=T)-theta[i]
    RMSE[i,] <- apply(re[re$select_index==i,name_point], 2, function(x) sqrt(mean((x-theta[i])^2, na.rm=T)))
    CP[i,] <- c(mean(apply(cbind(re[re$select_index==i,"stage2_L"] < theta[i] , re[re$select_index==i,"stage2_U"] > theta[i]),1,all),na.rm=T),
                mean(apply(cbind(re[re$select_index==i,"naive_L"] < theta[i] , re[re$select_index==i,"naive_U"] > theta[i]),1,all),na.rm=T),
                mean(apply(cbind(re[re$select_index==i,"RB-L"] < theta[i] , re[re$select_index==i,"RB-U"] > theta[i]),1,all),na.rm=T),
                mean(apply(cbind(re[re$select_index==i,"per-RB-L"] < theta[i] , re[re$select_index==i,"per-RB-U"] > theta[i]),1,all),na.rm=T))
    Width[i,] <- c(mean(re[re$select_index==i,"stage2_U"] - re[re$select_index==i,"stage2_L"] ,na.rm=T),
                   mean(re[re$select_index==i,"naive_U"] - re[re$select_index==i,"naive_L"] ,na.rm=T),
                   mean(re[re$select_index==i,"RB-U"] - re[re$select_index==i,"RB-L"] ,na.rm=T),
                   mean(re[re$select_index==i,"per-RB-U"] - re[re$select_index==i,"per-RB-L"] ,na.rm=T))
  }
  all_bias <- rbind(all_bias,bias)
  all_RMSE <- rbind(all_RMSE,RMSE)
  all_CP <- rbind(all_CP,CP)
  all_Width <- rbind(all_Width,Width)
}
write.csv(all_bias,paste("summary/",label,"-all_bias",".csv",sep=""))
write.csv(all_RMSE,paste("summary/",label,"-all_RMSE",".csv",sep=""))
write.csv(all_CP,paste("summary/",label,"-all_CP",".csv",sep=""))
write.csv(all_Width,paste("summary/",label,"-all_Width",".csv",sep=""))

##### plot #####
rm(list = ls())
# label <- "normal"
# label <- "binary"
# label <- "gamma"
for(label in c("normal","binary","gamma")){
  setwd("E:/0621/桌面/0919/new-res/2 stage/summary/")
  da_bias <- read.csv(paste(label,"-all_bias.csv",sep=""))
  da_RMSE <- read.csv(paste(label,"-all_RMSE.csv",sep=""))
  da_CP <- read.csv(paste(label,"-all_CP.csv",sep=""))
  da_Width <- read.csv(paste(label,"-all_Width.csv",sep=""))
  colnames(da_bias) <- colnames(da_RMSE) <- colnames(da_CP) <- colnames(da_Width) <- c("X","stage-wise","niave","RB","permutation")
  da_bias$select_index <- da_RMSE$select_index <- da_CP$select_index <- da_Width$select_index <- as.factor(rep(c(1:3),7))
  da_bias$power <- da_RMSE$power <- da_CP$power <- da_Width$power <- as.factor(rep(seq(0.3,0.9,0.1),each=3))
  
  ## plot
  da_bias$measure <- "conditional bias"
  da_RMSE$measure <- "conditional RMSE"
  da_CP$measure <- "conditional CP"
  da_Width$measure <- "conditional Width"
  da <- reshape2::melt(rbind(da_bias,da_RMSE,da_CP,da_Width),id.vars=c("select_index","measure","power"),
                       c("niave","stage-wise","RB","permutation"),
                       value.name="value",variable="method")
  da$measure <- factor(da$measure,levels = c("conditional bias","conditional RMSE","conditional CP","conditional Width"))
  da$treatment <- factor(da$select_index,c(1,2,3),c("Treatment 1","Treatment 2","Treatment 3"))
  ref <- data.frame(measure = "conditional CP", yintercept = 0.95)
  ref$measure <- factor(ref$measure,levels = c("conditional bias","conditional RMSE","conditional CP","conditional Width"))
  library(ggplot2)
  pdf(file=paste(label,"-result.pdf",sep=""),width=11,height=11)
  print(ggplot(da, aes(power,value,colour = method,group=method)) + 
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
          facet_grid(measure~treatment,scales="free_y"))
  dev.off()
  
  ## plot diff
  setwd("E:/0621/桌面/0919/new-res/2 stage/summary/")
  da_bias$diff <- da_bias$permutation-da_bias$RB
  da_RMSE$diff <- da_RMSE$permutation-da_RMSE$RB
  da_CP$diff <- da_CP$permutation-da_CP$RB
  da_Width$diff <- da_Width$permutation-da_Width$RB
  da <- rbind(da_bias,da_RMSE,da_CP,da_Width)
  da$measure <- factor(da$measure,levels = c("conditional bias","conditional RMSE","conditional CP","conditional Width"))
  da$treatment <- factor(da$select_index,c(1,2,3),c("Treatment 1","Treatment 2","Treatment 3"))
  ref <- data.frame(measure = c("conditional bias","conditional RMSE","conditional CP","conditional Width"), yintercept = 0)
  ref$measure <- factor(ref$measure,levels = c("conditional bias","conditional RMSE","conditional CP","conditional Width"))
  pdf(file=paste(label,"-result-diff.pdf",sep=""),width=11,height=11)
  print(ggplot(da, aes(power,diff,group = 1)) + 
          geom_point(size = 2) +   
          geom_line(linewidth = 1) +
          theme_bw()+
          theme(
            axis.title.x = element_text(size = 12, face = "bold"),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 13),
            legend.key.size = unit(0.8, "cm"),
            legend.position = "bottom",
            strip.text.x = element_text(size = 12, face = "bold"),
            strip.text.y = element_text(size = 12, face = "bold")
          )+
          ylab(NULL)+
          geom_hline(data = ref, aes(yintercept = yintercept), linetype = "dashed", color = "black")+
          facet_grid(measure~treatment,scales="free_y")
  )
  dev.off()
}
