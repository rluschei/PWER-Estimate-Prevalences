library(ggplot2)
library(tidyr)
library(MASS)
library(gridExtra)

### Export function

export <- function(m,N,i,df_temp){
  df_temp <- as.data.frame(df_temp)
  df_temp$rate <- ifelse(seq_len(nrow(df_temp)) %% 2 == 0, "fwer", "pwer")
  if(length(m)>1){
    names(df_temp) <- c(paste("m=", as.character(m), sep=""),"rate")
    boxplotsm(df_temp, i)
  }
  if(length(N)>1){
    names(df_temp) <- c(paste("N=", as.character(N), sep=""),"rate")
    boxplotsN(df_temp,i)
  }
  summary_table(df_temp[df_temp$rate=="pwer",-ncol(df_temp)], i, "pwer")
  summary_table(df_temp[df_temp$rate=="fwer",-ncol(df_temp)], i, "fwer")
}


## Boxplot functions
boxplotsm <- function(df_input, i){
  df_long <- gather(df_input, key="m", value="value", -rate)
  plot_pwer <- ggplot(df_long[df_long$rate=="pwer",], aes(x=m,y=value)) + geom_boxplot() + #coord_cartesian(ylim= c(0.0225,0.0275)) + 
    labs(x="number of biomarkers",y=expression(PWER(hat(c)))) +
    theme(text = element_text(size = 20), axis.title.x = element_text(margin=margin(t=10)))  
  ggsave(paste0("Plot",i,"pwer.pdf"), plot_pwer, height=8.27, width=11.69)
  
  plot_fwer <- ggplot(df_long[df_long$rate=="fwer",], aes(x=m,y=value)) + geom_boxplot() +
    labs(x="number of biomarkers",y=expression(max[J %subseteq% I]~FWER[J](hat(c)))) + 
    theme(text = element_text(size = 20), axis.title.x = element_text(margin=margin(t=10)))   
  ggsave(paste0("Plot",i,"fwer.pdf"), plot_fwer, height=8.27, width=11.69)
}


boxplotsN <- function(df_input, i){
  df_long <- gather(df_input, key="N", value="value", -rate)
  df_long$N <- factor(df_long$N, levels=c("N=25","N=50","N=100","N=150","N=200",
                                          "N=500"), ordered=T) # Set right order
  plot_pwer <- ggplot(df_long[df_long$rate=="pwer",], aes(x=N,y=value)) + geom_boxplot() +
    labs(x="overall sample size",y=expression(PWER(hat(c)))) + 
    theme(text = element_text(size = 20), axis.title.x = element_text(margin=margin(t=10)))
  ggsave(paste0("Plot",i,"pwer.pdf"), plot_pwer, height=8.27, width=11.69)
  
  plot_fwer <- ggplot(df_long[df_long$rate=="fwer",], aes(x=N,y=value)) + geom_boxplot() +
    labs(x="overall sample size",y=expression(max[J %subseteq% I]~FWER[J](hat(c)))) + 
    theme(text = element_text(size = 20), axis.title.x = element_text(margin=margin(t=10)))
  ggsave(paste0("Plot",i,"fwer.pdf"), plot_fwer, height=8.27, width=11.69)
}

boxplots_min <- function(df_input, i){
  df_input$group <- c("PWER","PWER_min","FWER","FWER_min")
  df_sub1 <- df_input[1:2,]
  df_long1 <- gather(df_sub1, key = "variable", value = "value", -group)
  df_sub2 <- df_input[3:4,]
  df_long2 <- gather(df_sub2, key = "variable", value = "value", -group)
  plot1 <- ggplot(df_long1, aes(x = group, y = value)) + 
    geom_boxplot() + xlab("") + ylab("") + 
    scale_x_discrete(labels=c(expression(PWER(hat(c))), expression(PWER(hat(c)[min])))) +
    theme(text = element_text(size = 20)) 
  plot2 <- ggplot(df_long2, aes(x = group, y = value)) +
    geom_boxplot() + xlab("") + ylab("") +
    scale_x_discrete(labels=c(expression(max[J %subseteq% I]~FWER[J](hat(c))), expression(max[J %subseteq% I]~FWER[J](hat(c)[min])))) +
    theme(text = element_text(size = 20)) 
  plot <- grid.arrange(plot1, plot2, ncol = 2)
  ggsave(paste0("Plot_min",i,".pdf"), plot, height=6, width=11.69)
}


## Summary table functions
summary_table <- function(df_input, i,value){
  summary_df <- data.frame(
    Mean = round(colMeans(df_input,na.rm=TRUE),digits=5),
    SD = round(apply(df_input, 2, sd,na.rm=TRUE),digits=5),
    Min = round(apply(df_input, 2, min,na.rm=TRUE),digits=5),
    Q1 = round(apply(df_input, 2, quantile,na.rm=TRUE, probs = 0.25),digits=5),
    Med = round(apply(df_input, 2, median,na.rm=TRUE),digits=5),
    Q3 = round(apply(df_input, 2, quantile, probs = 0.75,na.rm=TRUE),digits=5),
    Max = round(apply(df_input, 2, max,na.rm=TRUE),digits=5)
  )
  write.table(summary_df, file = paste0("summary",i,value,".txt"), sep="&", quote = F)    
}

summary_table_min <- function(df_input, i){
  summary_df <- data.frame(
    Mean = round(rowMeans(df_input,na.rm=TRUE),digits=5),
    SD = round(apply(df_input, 1, sd,na.rm=TRUE),digits=5),
    Min = round(apply(df_input, 1, min,na.rm=TRUE),digits=5),
    Q1 = round(apply(df_input, 1, quantile, probs = 0.25,na.rm=TRUE),digits=5),
    Med = round(apply(df_input, 1, median,na.rm=TRUE),digits=5),
    Q3 = round(apply(df_input, 1, quantile, probs = 0.75,na.rm=TRUE),digits=5),
    Max = round(apply(df_input, 1, max,na.rm=TRUE),digits=5)
  )
  write.table(summary_df, file = paste0("summary_min",i,".txt"), sep="&", quote = F) 
}
