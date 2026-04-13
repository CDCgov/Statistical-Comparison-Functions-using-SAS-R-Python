#-----------------------------Comparing 2 proportions by mid-p method------------------------------------------------------
# R script translated from SAS macro by Minn M. Soe, NHSN.
# Numerical accuracy of statistical results up to 6 decimal places was compared and validated with those from the corresponding SAS macro by Becky Lien, NHSN.
# Function: to analyse any simple 2 x 2 contingency table that shows the findings in two independent groups(within a single stratum). The data may be derived from an observational study (cross-sectional or cohort) or from a trial. 
# Use case scenarios: For hypothesis testing of comparing 2 proportions by Mid-p method based on hypergeometric distribution.


# ******************************************************************************

# Data input layout:
#   (E=exposure, D=disease)
#   cell A= E+,D+  
#   cell B= E-,D+
#   cell C= E+,D-
#   cell D= E-,D-
#                                     Disease/Outcome								   
#                            |    Yes                No       |                
#                 -----------|--------------------------------|-------         
#           Exposure 	  Yes	 |     A                  C			  |                
#                       No   |     B                  D			  |                
#                 -----------|--------------------------------|-------         
#           
# Output: mid-P value (2 sided)

# Note: No results will be generated if entry values are missing or impossible values. 

# ******************************************************************************
# Reference: 
# 1.David G. Kleinbaum, Lawrence L. Kupper. Epidemiologic Research: Principles and Quantitative Methods. 
# 2.Dean AG, Sullivan KM, Soe Minn Minn. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com, updated 2013/03/21
# 3.John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
# 4.Bernard Rosner. Fundamentals of Biostatistics' (5th edition) (Example 10.20, page 375). Two-tailed p-value calculated 
# as 2 times whichever is smallest: left-tail, right-tail, or 0.5. It tends to agree closely with Yates Chi-Square p-value.
# 5.David Martin, Harland Austin. An Efficient Program for Computing Conditional Maximum Likelihood Estimates and Exact Confidence Limits for a Common Odds Ratio. Epidemiology Vol. 2, No. 5 (Sep., 1991), pp. 359-362 

# ******************************************************************************

# Load necessary libraries
# library(tidyverse)
library(dplyr)
library(tidyr)

# Define the function as before
twobytwo <- function(A, B, C, D) {
  
  # Handle missing values
  if (all(c(A, B, C, D) == 0) || any(is.na(c(A, B, C, D))) || any(c(A, B, C, D) < 0)) {
    return(data.frame(Mid_P = NA))
  }
  
  # Input parameters
  Cell_A <- A
  Cell_B <- B
  Cell_C <- C
  Cell_D <- D
  Cell_r1 <- Cell_A + Cell_B
  Cell_r2 <- Cell_C + Cell_D
  Cell_c1 <- Cell_A + Cell_C
  Cell_c2 <- Cell_B + Cell_D
  t <- Cell_A + Cell_B + Cell_C + Cell_D
  
  # Constants
  ln_pi2 <- log(2 * pi)
  
  # Calculate LnProb1
  ln_fact <- function(z) {
    if (z < 2) return(0)
    if (z < 17) {
      f <- z
      while (z > 2) {
        z <- z - 1
        f <- f * z
      }
      return(log(f))
    } else {
      return((z + 0.5) * log(z) - z + ln_pi2 / 2 + 1 / (12 * z) - 1 / (360 * z^3) + 
               1 / (1260 * z^5) - 1 / (1680 * z^7))
    }
  }
  
  E1 <- ln_fact(Cell_r1)
  E2 <- ln_fact(Cell_r2)
  E3 <- ln_fact(Cell_c1)
  E4 <- ln_fact(Cell_c2)
  E5 <- ln_fact(t)
  LnProb1 <- E1 + E2 + E3 + E4 - E5
  
  # Initialize variables
  FisherP <- 0
  LeftP <- 0
  LeftP1 <- 0
  RightP <- 0
  RightP1 <- 0
  
  # Exact tests
  LoSlop <- min(Cell_A, Cell_D)
  HiSlop <- min(Cell_B, Cell_C)
  k <- Cell_A - LoSlop
  
  while (k <= Cell_A + HiSlop) {
    E1 <- ln_fact(k)
    E2 <- ln_fact(Cell_r1 - k)
    E3 <- ln_fact(Cell_c1 - k)
    E4 <- ln_fact(k + Cell_r2 - Cell_c1)
    P <- exp(LnProb1 - E1 - E2 - E3 - E4)
    
    if (k <= Cell_A) LeftP <- LeftP + P
    if (k < Cell_A) LeftP1 <- LeftP1 + P
    midp_left <- (LeftP - LeftP1) * 0.5 + LeftP1
    
    if (k >= Cell_A) RightP <- RightP + P
    if (k > Cell_A) RightP1 <- RightP1 + P
    midp_right <- (RightP - RightP1) * 0.5 + RightP1
    
    k <- k + 1
  }
  
  FisherP <- 2 * min(LeftP, RightP)
  Mid_P <- 2 * min(midp_left, midp_right)
  Mid_P <- max(0, min(1, Mid_P))  # Ensure Mid_P is between 0 and 1
  
  return(data.frame(Mid_P = Mid_P))
}



# --------------------------Example dataset----------------------------
df <- data.frame(
  id= c(1,   2, 3,  4),
  A = c(10, 20, 3, -3),
  B = c(20, 60, 5,  5),
  C = c(11, 80, 7, NA),
  D = c(21, 100, NA, 7)
)

# Apply the function row by row using dplyr
df1 <- df %>%
  rowwise() %>%
  mutate(result = list(twobytwo(A, B, C, D))) %>%
  unnest_wider(result)#%>% 
  #mutate(across(c(Mid_P), ~ round(., 3)))    #format to 3 decimal places

# Display the results
head(df1)
