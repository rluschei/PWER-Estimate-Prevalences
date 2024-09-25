
# Load libraries
source("libraries.R")

# Create a list of argument lists for the different cases to be treated
base_args <- list(case_no = 0, m = 2:8, N = 500, alpha = 0.025, indep_bio = TRUE, 
                  prev_mode = "random", distr = "t", strategy = "equal", 
                  single = FALSE, est = "MLE", het_var = F)

cases <- list(
  list(case_no = 1),
  list(case_no = 2, indep_bio = FALSE),
  list(case_no = 3, prev_mode = "fixed_equal"),
  list(case_no = 4, prev_mode = "fixed_0.5"),
  list(case_no = 5, est = "marg_est"),
  list(case_no = 6, distr = "z"),
  list(case_no = 7, single = TRUE),
  list(case_no = 8, strategy = "random"),
  list(case_no = 9, alpha = 0.01),
  list(case_no = 10, het_var = T),
  list(case_no = 11, m = 3, N = c(25, 50, 100, 150, 200, 500))
)

arg_lists <- lapply(cases, function(update) modifyList(base_args, update))

# Export list of cases
vec_as_char <- function(args) {  # function to handle vectors as single strings
  args$m <- paste(args$m, collapse = ", ")
  args$N <- paste(args$N, collapse = ", ")
  as.data.frame(t(unlist(args)), stringsAsFactors = FALSE)
}
arg_lists_df <- bind_rows(lapply(arg_lists, vec_as_char))
write.csv(arg_lists_df, file = "Cases_overview.csv", row.names = FALSE)