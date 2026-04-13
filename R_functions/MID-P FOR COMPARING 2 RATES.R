#-------------------mid-P method to compare 2 incidence rates---------------- 
# R script translated from SAS macro by Minn M. Soe, NHSN.
# Numerical accuracy of statistical results up to 6 decimal places was compared and validated with those from the corresponding SAS macro by Becky Lien, NHSN.
# Note: Outputs are mid-P value (two-sided) & optionally, point estimate and 95%CI of rate ratio.

# ******************************************************************************
# Data input parameters:
#   &O1= observed events in group-1
#   &PT1=person-time in group-1
#   &O2= observed events in group-2
#   &PT2= person-time in group-2
   
# OUTPUT:
#   mid-P=significance test for comparing between 2 rates
#   Rate Ratio= Rate in group-2 / Rate in group-1
#   LL=lower limit of 95%CI of rate ratio
#   UL=Upper limit of 95%CI of rate ratio

# Note: No results will be generated if entry values are missing, impossible values, or Person-time is not a whole number
# Output 'Rate Ratio' means RATE2/RATE1, ie, (O2/PT2)/(O1/PT1) assuming RATE2>RATE1. Therefore, if O1=0 then denominator is 0 and midP=NA

# ******************************************************************************
# Reference 
# 1. Erich L. Lehmann, Joseph P. Romano. Testing Statistical Hypotheses. Wiley, New York, 1959(Third edition, 2005, Springer Texts in Statistics) 
# 2. John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
# 3. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com, updated 2013/03/21
# 4. Austin, Harland (Emory). Epidemiology method-II: Statistical issues for Density type follow-up studies.

# ******************************************************************************
    
# Load necessary libraries
# library(tidyverse)
library(dplyr)
library(tidyr)

# Define functions equivalent to BinP and BinP2
binP <- function(x1, x2, PP, N) {
  q <- PP / (1 - PP)
  k <- 0
  v <- 1
  s <- 0
  tot <- 0
  while (k <= N) {
    tot <- tot + v
    if (k >= x1 && k <= x2) {s <- s + v}
    if (tot > 10^30) {
      s <- s / 10^30
      tot <- tot / 10^30
      v <- v / 10^30
    }
    k <- k + 1
    v <- v * q * (N + 1 - k) / k
  }
  return(s / tot)
}

binP2 <- function(x1, x2, PP, N) {
  q <- PP / (1 - PP)
  k <- 0
  v <- 1
  s <- 0
  tot <- 0
  while (k <= N) {
    tot <- tot + v
    if (k >= x1 && k <= x2) {s <- s + v}
    if (tot > 10^30) {
      s <- s / 10^30
      tot <- tot / 10^30
      v <- v / 10^30
    }
    k <- k + 1
    v <- v * q * (N + 1 - k) / k
  }
  return((s / tot) * 0.5)
}

# Define the main function for binomial calculations
  binom_rate <- function(O1, PT1, O2, PT2) {
    
  # Initialize result variables
  midP <- NA
  RATIO <- NA
  LL <- NA
  UL <- NA
  
  # Input validation
  #exclude if missing, impossible values or relations, PT is not a whole number
  if ((O1==0 && O2==0) || is.na(O1) || is.na(PT1) || is.na(O2) || is.na(PT2) || O1<0 || O2<0 || PT1<=0 || PT2<=0|| O1>PT1|| O2>PT2 || PT1 %% 1 != 0 || PT2 %% 1 != 0) {
    return(list(midP = midP, RATIO = RATIO, LL = LL, UL = UL))
  }
  
  vN <- O2 + O1
  vP <- O2 / vN
  
  # MID-P value for hypothesis testing
  RATIO1 <- O1 / PT1
  RATIO2 <- O2 / PT2
  O <- O1 + O2
  T <- PT1 + PT2
  
  # Function to calculate cumulative binomial probability
  prob_binomial <- function(p, n, x) {
    return(pbinom(x, n, p))
  }
  
  if (RATIO1 >= RATIO2) {
    p1 <- 1 - prob_binomial(PT1 / T, O, O1)
    p2 <- 0.5 * (prob_binomial(PT1 / T, O, O1) - prob_binomial(PT1 / T, O, O1 - 1))
    p3 <- p1 + p2
  } else {
    p1 <- 1 - prob_binomial(PT2 / T, O, O2)
    p2 <- 0.5 * (prob_binomial(PT2 / T, O, O2) - prob_binomial(PT2 / T, O, O2 - 1))
    p3 <- p1 + p2
  }
  
  midP <- 2 * p3
  midP <- min(max(midP, 0), 1)
  #midP <- ifelse(midP > 1, 1, ifelse(midP < 0, 0, midP))
  
  # RATIO and Interval Estimation (if O1=0, no ratio & CI computed)
  if (O1 != 0) {
    RATIO <- (O2 / PT2) / (O1 / PT1)
    
    # Lower Limit
    if (O2 == 0) {
      DL <- 0
    } else {
      p2 <- vP / 2
      vsL <- 0
      vsH <- vP
      p <- 2.5 / 100
      while ((vsH - vsL) > 10^(-5)) {
        BinP <- binP(x1 = O2 + 1, x2 = vN, PP = p2, N = vN)
        BinP2 <- binP2(x1 = O2, x2 = O2, PP = p2, N = vN)
        if ((BinP + BinP2) > p) {
          vsH <- p2
          p2 <- (vsL + p2) / 2
        } else {
          vsL <- p2
          p2 <- (vsH + p2) / 2
        }
      }
      DL <- p2
    }
    
    # Upper Limit
    if (O2 == vN) {
      DU <- 1
    } else {
      p3 <- (1 + vP) / 2
      vsL <- vP
      vsH <- 1
      p <- 2.5 / 100
      while ((vsH - vsL) > 10^(-5)) {
        BinP <- binP(x1 = 0, x2 = O2 - 1, PP = p3, N = vN)
        BinP2 <- binP2(x1 = O2, x2 = O2, PP = p3, N = vN)
        if ((BinP + BinP2) < p) {
          vsH <- p3
          p3 <- (vsL + p3) / 2
        } else {
          vsL <- p3
          p3 <- (vsH + p3) / 2
        }
      }
      DU <- p3
    }
    
    LL <- (DL * PT1) / ((1 - DL) * PT2)
    UL <- (DU * PT1) / ((1 - DU) * PT2)

    if (RATIO == 0) {LL <- NA}  #set LL to missing if RATIO=0
  } 
  
  else {
    RATIO <- NA
    LL <- NA
    UL <- NA
  }
  
  return(data.frame(midP = midP, RATIO = RATIO, LL = LL, UL = UL))
}

  
# --------------------------Example dataset---------------------------------------------
#comparing 2 rates between 2 groups
# obs=number of events,  pt=person-time
    
df <- data.frame(
  id      = c(1,2,3,4),
  obs2017 = c(NA,1,1,1),
  pt2017  = c(2.3,2,2,2),
  obs2016 = c(2,3,3,0),
  pt2016  = c(4.5,-4,4,5)
)

# Apply the binom function row-wise
df1 <- df %>%
  rowwise() %>%
  mutate(binom_result = list(binom_rate(O1 = obs2017, PT1 = pt2017, O2 = obs2016, PT2 = pt2016))) %>%
  unnest_wider(binom_result)#%>% 
  # mutate(across(c(midP, RATIO, LL, UL), ~ round(., 3)))    #format to 3 decimal places

# Display the results
head(df1)


