
# Load libraries
source("libraries.R")

# Load functions and list of cases
source("functions.R")
source("cases.R")

# Boxplot function
boxplot <- function(df, case_no){
  if(case_no < 10) {
    key <- "m"
    df_long <- gather(df, key = key, value = "value")
    x_label <- "number of biomarkers m"
  }
  if(case_no == 10) {
    key <- "N"
    df_long <- gather(df, key = key, value = "value")
    df_long$key <- factor(df_long$key, levels = c("25","50","100","150","200","500"), ordered = TRUE)
    x_label <- "sample size N"
  }
  y_label <- expression(PWER(hat(c)))
  
  plot <- ggplot(df_long, aes(x = key, y = value)) + geom_boxplot() +
    labs(x = x_label, y = y_label) + 
    theme(text = element_text(size = 20), axis.title.x = element_text(margin = margin(t = 10)))
  ggsave(paste0("Plot_pwer_case", case_no, ".pdf"), plot, height = 8.27, 
         width = 11.69, path = "Results")
}

# Wrapper function that creates data frames with the simulated values and calls 
# the boxplot and summary table function
wrap <- function(args) {
  df <- as.data.frame(mapply(sim, m = args$m, N = args$N, alpha = args$alpha, 
                             indep_bio = args$indep_bio, prev_mode = args$prev_mode, 
                             distr = args$distr, strategy = args$strategy, 
                             tmt_comp = args$tmt_comp, est = args$est, error_rate = "pwer"))
  names(df) <- ifelse(rep(args$case_no, length.out = ncol(df)) < 10, args$m, args$N)
  assign(paste0("df_pwer", args$case_no), df, envir = .GlobalEnv)
  boxplot(df, args$case_no)
  summary_tbl(df, "pwer", args$case_no)
}

# Run simulations 
invisible(map(arg_lists, function(args) wrap(args)))









