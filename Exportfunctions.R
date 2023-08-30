## Boxplot functions
boxplots <- function(df_input, i, y, value){
  df_long <- gather(df_input, key="m", value="pwer")
  if(value=="pwer"){
    plot <- ggplot(df_long, aes(x=m,y=pwer)) + geom_boxplot(outlier.shape=NA) +
      labs(x="",y=expression(PWER(hat(c)))) + coord_cartesian(ylim=y) +
      theme(text = element_text(size = 20))   
  }
  if(value=="fwer"){
    plot <- ggplot(df_long, aes(x=m,y=pwer)) + geom_boxplot(outlier.shape=NA) +
      labs(x="",y=expression(max[J %subseteq% I]~FWER[J](hat(c)))) + 
      coord_cartesian(ylim=y) + theme(text = element_text(size = 20))   
  }
  ggsave(paste0("Plot",i,".pdf"), plot, height=8.27, width=11.69)
}

boxplotsN <- function(df_input, i, y, value){
  df_long <- gather(df_input, key="N", value="pwer")
  df_long$N <- factor(df_long$N, levels=c("N=25","N=50","N=100","N=150","N=200",
                                          "N=500"), ordered=T) # Set right order
  if(value=="pwer"){
    plot <- ggplot(df_long, aes(x=N,y=pwer)) + geom_boxplot(outlier.shape=NA) +
      labs(x="",y=expression(PWER(hat(c)))) + coord_cartesian(ylim=y) +
      theme(text = element_text(size = 20))
  }
  if(value=="fwer"){
    plot <- ggplot(df_long, aes(x=N,y=pwer)) + geom_boxplot(outlier.shape=NA) +
      labs(x="",y=expression(max[J %subseteq% I]~FWER[J](hat(c)))) + 
      coord_cartesian(ylim=y) + theme(text = element_text(size = 20))
  }
  ggsave(paste0("Plot",i,".pdf"), plot, height=8.27, width=11.69)
}

boxplots_min <- function(df_input, i){
  df_input$group <- c("PWER","PWER_min","FWER","FWER_min")
  df_sub1 <- df_input[1:2,]
  df_long1 <- gather(df_sub1, key = "variable", value = "value", -group)
  df_sub2 <- df_input[3:4,]
  df_long2 <- gather(df_sub2, key = "variable", value = "value", -group)
  plot1 <- ggplot(df_long1, aes(x = group, y = value)) + 
    geom_boxplot(outlier.shape=NA) + xlab("") + ylab("") + 
    scale_x_discrete(labels=c(expression(PWER(hat(c))), expression(PWER(hat(c)[min])))) +
    theme(text = element_text(size = 20)) 
  plot2 <- ggplot(df_long2, aes(x = group, y = value)) +
    geom_boxplot(outlier.shape=NA) + xlab("") + ylab("") +
    scale_x_discrete(labels=c(expression(max[J %subseteq% I]~FWER[J](hat(c))), expression(max[J %subseteq% I]~FWER[J](hat(c)[min])))) +
    theme(text = element_text(size = 20)) 
  plot <- grid.arrange(plot1, plot2, ncol = 2)
  ggsave(paste0("Plot",i,".pdf"), plot, height=6, width=11.69)
}


## Summary table functions
summary_table <- function(df_input, i){
  summary_df <- data.frame(
    Mean = round(colMeans(df_input),digits=5),
    SD = round(apply(df_input, 2, sd),digits=5),
    Min = round(apply(df_input, 2, min),digits=5),
    Q1 = round(apply(df_input, 2, quantile, probs = 0.25),digits=5),
    Med = round(apply(df_input, 2, median),digits=5),
    Q3 = round(apply(df_input, 2, quantile, probs = 0.75),digits=5),
    Max = round(apply(df_input, 2, max),digits=5)
  )
  write.table(summary_df, file = paste0("summary",i,".txt"), sep="&", quote = F)    
}

summary_table_min <- function(df_input, i, Nsim){
  summary_df <- data.frame(
    N = Nsim-rowSums(is.na(df_input)),
    Mean = round(rowMeans(df_input,na.rm=TRUE),digits=5),
    SD = round(apply(df_input, 1, sd,na.rm=TRUE),digits=5),
    Min = round(apply(df_input, 1, min,na.rm=TRUE),digits=5),
    Q1 = round(apply(df_input, 1, quantile, probs = 0.25,na.rm=TRUE),digits=5),
    Med = round(apply(df_input, 1, median,na.rm=TRUE),digits=5),
    Q3 = round(apply(df_input, 1, quantile, probs = 0.75,na.rm=TRUE),digits=5),
    Max = round(apply(df_input, 1, max,na.rm=TRUE),digits=5)
  )
  write.table(summary_df, file = paste0("summary",i,".txt"), sep="&", quote = F) 
}
