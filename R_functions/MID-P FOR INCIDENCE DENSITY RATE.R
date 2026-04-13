#------------------mid-P method to compute 95%CI of single incidence rate-------------------- 
# R script translated from SAS macro by Minn M. Soe, NHSN.
# Numerical accuracy of statistical results up to 6 decimal places was compared and validated with those from the corresponding SAS macro by Becky Lien, NHSN.

# ****************************************************************************** 
# Input parameters:  numer=numerator,  denom=denominator.

# output: incidence density rate per 1000 & lower and upper limit of 95%CI (rate_l & rate_u).

# Note: No results will be generated if entry values are missing or impossible values. 

# ****************************************************************************** 
# Reference 
# 1. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version 2.3.1. www.OpenEpi.com
# 2. Rothman KJ,  Boice JD  Jr:  Epidemiologic analysis with a programmable calculator.   NIH Pub No. 79-1649.  Bethesda, MD:  National Institutes of Health, 1979;31-32.
# 3. Rothman KJ, Greenland S.  Modern Epidemiology, 2nd Edition.  Lippincott-Raven Publishers, Philadelphia, 1998.

# ****************************************************************************** 

# Load necessary libraries
# library(tidyverse)
library(dplyr)
library(tidyr)

# Fish function
fish <- function(z, x1, x2) {
  q <- 1
  tot <- 0
  s <- 0
  k <- 0
  while (k < z || q > (tot * 10^(-10))) {
    tot <- tot + q
    if (k >= x1 && k <= x2) {s <- s + q}
    if (tot > 10^30) {
      s <- s / 10^30
      tot <- tot / 10^30
      q <- q / 10^30
    }
    k <- k + 1
    q <- q * z / k
  }
  return(s / tot)
}

# Fish2 function
fish2 <- function(z, x1, x2) {
  q <- 1
  tot <- 0
  s <- 0
  k <- 0
  while (k < z || q > (tot * 10^(-10))) {
    tot <- tot + q
    if (k >= x1 && k <= x2) {s <- s + q}
    if (tot > 10^30) {
      s <- s / 10^30
      tot <- tot / 10^30
      q <- q / 10^30
    }
    k <- k + 1
    q <- q * z / k
  }
  return(0.5 * (s / tot))
}

rateCIComp <- function(numer, denom) {
  rate<- (numer/denom)*1000
  
  # Check for NA values
  if (is.na(numer) || is.na(denom) || numer<0 || denom<=0 || denom<numer || denom %% 1 != 0) {
    return(list(rate=NA, rate_l = NA, rate_u = NA))
  }
  
  # Lower tail
  if (numer == 0) {
    rate_l <- NA
  } else {
    v <- 0.5
    dv <- 0.5
    p <- 2.5 / 100
    while (dv > 10^(-5)) {
      dv <- dv / 2
      poisP <- fish(z = (1 + numer) * v / (1 - v), x1 = (numer + 1), x2 = 10^10)
      poisP2 <- fish2(z = (1 + numer) * v / (1 - v), x1 = numer, x2 = numer)
      if ((poisP + poisP2) > p) {
        v <- v - dv
      } else {
        v <- v + dv
      }
    }
    rate_l <- (((1 + numer) * v / (1 - v)) / denom) * 1000
  }
    
  # Upper tail
  v <- 0.5
  dv <- 0.5
  p <- 2.5 / 100
  while (dv > 10^(-5)) {
    dv <- dv / 2
    poisP <- fish(z = (1 + numer) * v / (1 - v), x1 = 0, x2 = (numer - 1))
    poisP2 <- fish2(z = (1 + numer) * v / (1 - v), x1 = numer, x2 = numer)
    if ((poisP + poisP2) < p) {
      v <- v - dv
    } else {
      v <- v + dv
    }
  }
  rate_u <- (((1 + numer) * v / (1 - v)) / denom) * 1000
 
  return(list(rate=rate, rate_l = rate_l, rate_u = rate_u))
  }


# CREATE AN EXAMPLE DATAFRAME
df <- data.frame(
    id  = c(1, 2),
    obs = c(0,10),
    pt  = c(10,10)
)


# Apply the function row-wise and add multiple columns
df1 <- df %>%
  rowwise() %>%
  mutate(result = list(rateCIComp(obs,pt)))%>%  #replace with num&denom in dframe
  unnest_wider(result) 
# df1 <- df1 %>% mutate(across(c(rate, rate_l, rate_u), ~ format(., nsmall = 3)))#keep trailing zeros
# df1<-df1 %>% mutate(across(c(rate, rate_l, rate_u), ~ round(., 3)))    #format to 3 decimal places
# Display the results
head(df1)





