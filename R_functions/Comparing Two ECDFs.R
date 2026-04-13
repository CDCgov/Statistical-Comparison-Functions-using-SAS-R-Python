#------------------------------------------------------------------------------#
#   Comparing Two ECDFs.R                                                      #
#                                                                              #
#   Function:  compare_two_ecdfs                                               #
#                                                                              #
#   Inputs:   1) Grouping variable called group_var                            #
#             2) Group 1 label "group 1"                                       #
#             3) Group 2 label "group 2"                                       #
#                                                                              #
#   Outputs:   1) Brown-Mood Median and ECDF (Anderson-Darling, DTS) Tests     #
#              2) Basic and Advanced Smoothed ECDF Plots                       #
#                                                                              #
#  # Example usage:                                                            #
#   compare_two_ecdfs(                                                         #
#    data = df,                                                                #
#    outcome = "some_other_outcome",                                           #
#    group_var = "arm",                                                        #
#    group1 = "group_1",                                                       #
#    group2 = "group_2"                                                        #
#   )                                                                          #
#                                                                              #
#                                                                              #
#   Four examples have been included that display use of this function         #
#                                                                              #
#   Developed by Jonathan R. Edwards   Feb 2, 2026                             #
#------------------------------------------------------------------------------#


compare_two_ecdfs <- function(
    data,
    outcome,
    group_var = "surface",
    group1,
    group2,
    plot_label = "ECDF Comparison",
    min_n = 5
) {
  # Packages
  if (!requireNamespace("PMCMRplus", quietly = TRUE)) stop("Install PMCMRplus.")
  if (!requireNamespace("twosamples", quietly = TRUE)) stop("Install twosamples.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2.")
  
  library(PMCMRplus)
  library(twosamples)
  library(ggplot2)
  
  # Allow outcome to be passed as unquoted or quoted
  outcome <- if (is.character(outcome)) outcome else deparse(substitute(outcome))
  
  # Basic checks
  if (!group_var %in% names(data)) stop(sprintf("'%s' not found in data.", group_var))
  if (!outcome %in% names(data)) stop(sprintf("Outcome '%s' not found in data.", outcome))
  
  # Ensure grouping is a factor
  data[[group_var]] <- as.factor(data[[group_var]])
  
  # Extract groups, drop NAs
  g1 <- na.omit(data[[outcome]][data[[group_var]] == group1])
  g2 <- na.omit(data[[outcome]][data[[group_var]] == group2])
  
  n1 <- length(g1)
  n2 <- length(g2)
  
  # Validate input
  if (n1 < min_n || n2 < min_n || var(g1) == 0 || var(g2) == 0) {
    warning("Insufficient data or zero variance in one of the groups.")
    return(NULL)
  }
  
  # Tests
  bm_test <- medianTest(stats::as.formula(paste(outcome, "~", group_var)), data = data)
  ad_twosamples <- ad_test(g1, g2)
  
  dts_p <- tryCatch({
    dts <- dts_test(g1, g2)
    as.numeric(dts["P-Value"])
  }, error = function(e) NA_real_)
  
  bm_p <- bm_test$p.value
  ad_p <- as.numeric(ad_twosamples["P-Value"])
  
  test_summary <- paste(
    sprintf("n1=%d, n2=%d", n1, n2),
    sprintf("Brown-Mood Median p=%.3g", bm_p),
    sprintf("Anderson-Darling ECDF p=%.3g", ad_p),
    sprintf("DTS ECDF p=%.3g", dts_p),
    sep = "  |  "
  )
  

  # Plot data
  df <- data.frame(
    value = c(g1, g2),
    group = factor(c(rep(group1, n1), rep(group2, n2)),
                   levels = c(group1, group2))
  )
  
  surface_colors <- stats::setNames(
    c("royalblue3", "darkseagreen3"),
    c(group1, group2)
  )
  
  # Basic ECDF Plot
  p1 <- ggplot(df, aes(x = value, color = group)) +
    stat_ecdf(linewidth = 1.3) +
    labs(
      title = paste0("Empirical CDFs of ", outcome, ": ", plot_label),
      subtitle = test_summary,
      x = outcome, y = "ECDF", color = group_var
    ) +
    scale_color_manual(values = surface_colors) +
    theme_light(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
  
  print(p1)
  
  # ECDF Plot with Shaded Difference
  x_vals <- sort(unique(c(g1, g2)))
  ecdf_g1 <- ecdf(g1)
  ecdf_g2 <- ecdf(g2)
  
  df_ecdf <- data.frame(
    x = x_vals,
    ecdf_g1 = ecdf_g1(x_vals),
    ecdf_g2 = ecdf_g2(x_vals),
    ymin = pmin(ecdf_g1(x_vals), ecdf_g2(x_vals)),
    ymax = pmax(ecdf_g1(x_vals), ecdf_g2(x_vals))
  )
  
  p2 <- ggplot(df_ecdf, aes(x = x)) +
    geom_line(aes(y = ecdf_g1, color = group1), linewidth = 1.3) +
    geom_line(aes(y = ecdf_g2, color = group2), linewidth = 1.3) +
    geom_ribbon(aes(ymin = ymin, ymax = ymax),
                fill = "slategray4", alpha = 0.3) +
    labs(
      title = paste0("ECDFs of ", outcome, " with Shaded Difference: ", plot_label),
      subtitle = test_summary,
      x = outcome, y = "Empirical CDF", color = group_var
    ) +
    scale_color_manual(values = surface_colors) +
    theme_light(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    )
  
  print(p2)
  
  invisible(list(
    n1 = n1,
    n2 = n2,
    p_brown_mood = bm_p,
    p_anderson_darling = ad_p,
    p_dts = dts_p,
    subtitle = test_summary,
    plot_ecdf = p1,
    plot_ecdf_shaded = p2
  ))
}




#------------------------------------------------------------------------------#
#   Example:  Generic data frame from a normal distribution                    #
#             function: compare_two_ecdfs                                      #
#                                                                              #
#------------------------------------------------------------------------------#

set.seed(123)

example_df <- data.frame(
  subject_id = 1:120,
  group = rep(c("Group_1", "Group_2"), each = 60),
  outcome_value = c(
    rnorm(60, mean = 2.2, sd = 0.45),   # Group 1 distribution
    rnorm(60, mean = 2.7, sd = 0.50)    # Group 2 distribution
  )
)

# Introduce a few missing values (optional, realistic)
example_df$outcome_value[sample(1:120, 6)] <- NA

str(example_df)


compare_two_ecdfs(
  data = example_df,
  outcome = "outcome_value",
  group_var = "group",
  group1 = "Group_1",
  group2 = "Group_2",
  plot_label = "Example normal outcome"
)

