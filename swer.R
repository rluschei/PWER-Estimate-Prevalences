
# Load libraries
source("libraries.R")

# Load functions and list of cases
source("functions.R")
source("cases.R")

# Boxplot function
boxplot <- function(df_max_swer, df_mean_swer, case_no){
  if(case_no < 11) {
    key <- "m"
    df_max_swer_long <- gather(df_max_swer, key = key, value = "value")
    df_mean_swer_long <- gather(df_mean_swer, key = key, value = "value")
    x_label <- "number of biomarkers m"
  }
  if(case_no == 11) {
    key <- "N"
    df_max_swer_long <- gather(df_max_swer, key = key, value = "value")
    df_mean_swer_long <- gather(df_mean_swer, key = key, value = "value")
    df_max_swer_long$key <- factor(df_max_swer_long$key, levels = c("25","50","100","150","200","500"), ordered = TRUE)
    df_mean_swer_long$key <- factor(df_mean_swer_long$key, levels = c("25","50","100","150","200","500"), ordered = TRUE)
    x_label <- "sample size N"
  }
  
  df_max_swer_long$error_rate <- "maximal SWER"
  df_mean_swer_long$error_rate <- "mean SWER"
  
  combined_df <- rbind(df_max_swer_long, df_mean_swer_long)
  
  plot <- ggplot(combined_df, aes(x = key, y = value, fill = error_rate)) + 
    geom_boxplot(position = position_dodge(width = 0.75)) + 
    labs(x = x_label, y = "") + 
    scale_fill_manual(values = c("maximal SWER" = "skyblue", "mean SWER" = "lightgreen"), name = NULL) +
    theme(text = element_text(size = 20), axis.title.x = element_text(margin = margin(t = 10)))
  
  ggsave(paste0("Plot_swer_case", case_no, ".pdf"), plot, height = 8.27, width = 15, path = "Results")
}


# Wrapper function that creates data frames with the simulated values and calls 
# the boxplot and summary table function
wrap <- function(args) {
  dfs <- lapply(c("max_swer", "mean_swer"), function(error_rate) {
    df <- as.data.frame(mapply(sim, m = args$m, N = args$N, alpha = args$alpha, 
                               indep_bio = args$indep_bio, prev_mode = args$prev_mode, 
                               distr = args$distr, strategy = args$strategy, 
                               single = args$single, est = args$est, 
                               het_var = args$het_var, error_rate = error_rate))
    names(df) <- ifelse(rep(args$case_no, length.out = ncol(df)) < 11, args$m, args$N)
    assign(paste0("df_", error_rate, args$case_no), df, envir = .GlobalEnv)
  })
  
  boxplot(dfs[[1]], dfs[[2]], args$case_no)
  summary_tbl(dfs[[1]], "max_swer", args$case_no)
  summary_tbl(dfs[[2]], "mean_swer", args$case_no)
}

invisible(map(arg_lists, function(args) wrap(args)))









