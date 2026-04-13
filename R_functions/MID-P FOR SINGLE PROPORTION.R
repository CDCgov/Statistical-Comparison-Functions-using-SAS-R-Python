#----------------------mid-P method to compute 95%CI of a proportion-------------------- 
# R script translated from SAS macro by Minn M. Soe, NHSN.
# Numerical accuracy of statistical results up to 6 decimal places was compared and validated with those from the corresponding SAS macro by Becky Lien, NHSN.

# ********************************************************************************

# Input parameters:  vx=numerator,  vn=denominator.
 
# Output: vp=proportion, dl=lower limit of 95%CI, du=upper limit of 95%CI

# Note: No results will be generated if entry values are missing or impossible values. 
# The BinP and BinP2 are numerical methods for computing mid-P values via the likelihood ratio.

# ********************************************************************************
# Reference 
# 1. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version 2.3.1. www.OpenEpi.com
# 2. Rothman KJ,  Boice JD  Jr:  Epidemiologic analysis with a programmable calculator. NIH Pub No. 79-1649.  Bethesda, MD:  National Institutes of Health, 1979;31-32.
# 3. Fleiss JL.  Statistical Methods for Rates and Proportions, 2nd Ed. John Wiley & Sons, New York, 1981.
# 4. Newcombe RG.  Two-sided confidence intervals for the single proportion: comparison of seven methods.  Statistics in Medicine 1988;17:857-872.
# 5. Stein Vollset. CONFIDENCE INTERVALS FOR A BINOMIAL PROPORTION STATISTICS IN MEDICINE,12,809-824 (1993)*/
  
# ********************************************************************************
    
# Load necessary libraries
# library(tidyverse)
library(dplyr)
library(tidyr)

BinP <- function(x1, x2, p2, vN) {
  q <- p2 / (1 - p2)
  k <- 0
  v <- 1
  s <- 0
  tot <- 0
  while (k <= vN) {
    tot <- tot + v
    if (k >= x1 && k <= x2) s <- s + v
    if (tot > 1e30) {
      s <- s / 1e30
      tot <- tot / 1e30
      v <- v / 1e30
    }
    k <- k + 1
    v <- v * q * (vN + 1 - k) / k
  }
  return(s / tot)
}

BinP2 <- function(x1, x2, p2, vN) {
  q <- p2 / (1 - p2)
  k <- 0
  v <- 1
  s <- 0
  tot <- 0
  while (k <= vN) {
    tot <- tot + v
    if (k >= x1 && k <= x2) s <- s + v
    if (tot > 1e30) {
      s <- s / 1e30
      tot <- tot / 1e30
      v <- v / 1e30
    }
    k <- k + 1
    v <- v * q * (vN + 1 - k) / k
  }
  return((s / tot) * 0.5)
}

binom_midp_ci <- function(vX, vN) {
  # Initialize result values
  vp <- NA
  dl <- NA
  du <- NA
  
  # Handle missing/invalid input
  if (is.na(vX) || is.na(vN) || vX < 0 || vN <= 0 || vX > vN || vX %% 1 != 0 || vN %% 1 != 0) {
    return(list(vp = vp, dl = dl, du = du))
  }
  
  vp <- vX / vN
  
  ### Lower Limit ###
  if (vX == 0) {
    dl <- NA
  } else {
    p2 <- vp / 2
    vsL <- 0
    vsH <- vp
    p <- 0.025
    
 
    
    while ((vsH - vsL) > 1e-5) {
      b1 <- BinP(vX + 1, vN, p2, vN)
      b2 <- BinP2(vX, vX, p2, vN)
      if ((b1 + b2) > p) {
        vsH <- p2
        p2 <- (vsL + p2) / 2
      } else {
        vsL <- p2
        p2 <- (vsH + p2) / 2
      }
    }
    dl <- p2
  }
  
  ### Upper Limit ###
  if (vX == vN) {
    du <- 1
  } else {
    p2 <- (1 + vp) / 2
    vsL <- vp
    vsH <- 1
    p <- 0.025
    
    while ((vsH - vsL) > 1e-5) {
      b1 <- BinP(0, vX - 1, p2, vN)
      b2 <- BinP2(vX, vX, p2, vN)
      if ((b1 + b2) < p) {
        vsH <- p2
        p2 <- (vsL + p2) / 2
      } else {
        vsL <- p2
        p2 <- (vsH + p2) / 2
      }
    }
    du <- p2
  }
  
  return(list(vp = vp, dl = dl, du = du))
}


# CREATE AN EXAMPLE DATAFRAME
df <- data.frame(
  id   = c(1, 2, 3),
  obs  = c(1, 10, NA),
  total= c(2, 4, 0.03)
)

# Apply the function row-wise and add multiple columns
df1 <- df %>%
  rowwise() %>%
  mutate(result = list(binom_midp_ci(obs,total)))%>%  #replace with num&denom in dframe
  unnest_wider(result) #%>% 
  # mutate(across(c(vp, dl, du), ~ round(., 3)))    #format to 3 decimal places

# Display the results
head(df1)


