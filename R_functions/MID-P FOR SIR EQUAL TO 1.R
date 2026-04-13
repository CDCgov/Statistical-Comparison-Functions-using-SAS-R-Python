#----------------------------computing mid-p value & 95%CI of single SIR--------------------------
# R script translated from SAS macro by Minn M. Soe, NHSN.
# Numerical accuracy of statistical results up to 6 decimal places was compared and validated with those from the corresponding SAS macro by Lindsay Dunham-Maher, NHSN.

# ******************************************************************************  
# Input parameters:  OBS=observed value,  EXP=predicted value

# Output: SIR=SIR, LLIMIT & ULIMIT=lower and upper limit of 95%CI of SIR

# Note: No results will be generated if entry values are missing or impossible values.
# When SIR=0, the lower limit (LLIMIT) of 95%CI is set to missing by default. A user may change it to LLIMIT=0 if desired.

# ******************************************************************************  
# REFERENCE: 
# 1. JH ABRAMSON. WINPEPI PROGRAMS. DESCRIBE MANUAL (VERSION 2.42), PAGE-52. Available at 'http://www.brixtonhealth.com/pepi4windows.html'
# 2. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com.
# 3. Rothman KJ, Boice JD Jr: Epidemiologic analysis with a programmable calculator. NIH Pub No. 79-1649. Bethesda, MD: National Institutes of Health, 1979;31-32.
# 4. Rothman KJ, Greenland S.  Modern Epidemiology, 2nd Edition.  Lippincott-Raven Publishers, Philadelphia, 1998.
# 5. GEOFFREY RC, SHU-YING Y. MID-P CONFIDENCE INTERVALS FOR THE POISSON EXPECTATION. STATISTICS IN MEDICINE, 13,2189-2203 (1994)

# ******************************************************************************   

# Load necessary libraries
# library(tidyverse)
library(dplyr)
library(tidyr)

sir_function <- function(obs, pred) {
  e <- exp(1)  # Euler's number
  
  # Initialize results
  midp <- NA
  Byar_p <- NA
  sir_val <- NA
  ulimit <- NA
  llimit <- NA
  
  if (is.na(obs) || is.na(pred) || pred<1 || obs<0 || obs%%1 != 0) {
    return(list(midp = midp, sir = sir_val, llimit = llimit, ulimit = ulimit))
  }
  
  sir_val <- obs / pred
  
  if (obs <= 100) {
    total <- 0
    k <- obs - 1
    while (k >= 0) {
      numerator <- (e^(-pred)) * (pred^k)
      denom <- factorial(k)
      subtotal <- numerator / denom
      total <- total + subtotal
      k <- k - 1
    }
    
    num <- (e^(-pred)) * (pred^obs)
    deno <- factorial(obs)
    aa <- (num / deno) * 0.5
    
    if (obs < pred) {
      midp <- 2 * (aa + total)
    } else {
      midp <- 2 * (1 - (aa + total))
    }
    
    midp <- min(max(midp, 0), 1)
  }
  
  # When obs > 100, use Byar Poisson approximation
  if (obs > 100) {
    obs_ <- if (obs <= pred) obs + 1 else obs
    z <- sqrt(9 * obs_) * (1 - 1 / (9 * obs_) - (pred / obs_)^(1 / 3))
    
    Byar_p <- if (obs <= pred) {
      2 * pnorm(z)
    } else {
      2 * (1 - pnorm(z))
    }
    
    midp <- min(max(Byar_p, 0), 1)
  }
  
  # Functions for exact CI calculation
  fish <- function(z, x1, x2) {
    q <- 1; tot <- 0; s <- 0; k <- 0
    while (k < z || q > (tot * 10^(-10))) {
      tot <- tot + q
      if (k >= x1 && k <= x2) s <- s + q
      if (tot > 1e30) {
        s <- s / 1e30; tot <- tot / 1e30; q <- q / 1e30
      }
      k <- k + 1
      q <- q * z / k
    }
    return(s / tot)
  }
  
  fish2 <- function(z, x1, x2) {
    q <- 1; tot <- 0; s <- 0; k <- 0
    while (k < z || q > (tot * 10^(-10))) {
      tot <- tot + q
      if (k >= x1 && k <= x2) s <- s + q
      if (tot > 1e30) {
        s <- s / 1e30; tot <- tot / 1e30; q <- q / 1e30
      }
      k <- k + 1
      q <- q * z / k
    }
    return(0.5 * (s / tot))
  }
  
  # Lower confidence limit
  if (obs != 0) {
    v <- 0.5; dv <- 0.5; p <- 2.5 / 100
    while (dv > 1e-5) {
      dv <- dv / 2
      zval <- (1 + obs) * v / (1 - v)
      component_a <- fish(z = zval, x1 = (obs + 1), x2 = 1e10)
      component_b <- fish2(z = zval, x1 = obs, x2 = obs)
      if ((component_a + component_b) > p) {
        v <- v - dv
      } else {
        v <- v + dv
      }
    }
    llimit <- ((1 + obs) * v / (1 - v)) / pred
  }
  
  # Upper confidence limit
  v <- 0.5; dv <- 0.5; p <- 2.5 / 100
  while (dv > 1e-5) {
    dv <- dv / 2
    zval <- (1 + obs) * v / (1 - v)
    component_a <- fish(z = zval, x1 = 0, x2 = (obs - 1))
    component_b <- fish2(z = zval, x1 = obs, x2 = obs)
    if ((component_a + component_b) < p) {
      v <- v - dv
    } else {
      v <- v + dv
    }
  }
  ulimit <- ((1 + obs) * v / (1 - v)) / pred
  
  return(list(midp = midp, sir = sir_val, llimit = llimit, ulimit = ulimit))
}



#---------------------------EXAMPLE DATAFRAME-------------------------------------
df <- data.frame(
  id    = c(1, 2, 3),
  num   = c(1, 2, NA),
  denom = c(2, 0.6, 3.3)
)

#substitute with real data set & apply corresponding variable names for numerator & denominator in sir_function

# Apply the sir_function row-wise 
df1 <- df %>%
  rowwise() %>%
  mutate(result = list(sir_function(obs=num, pred=denom)))%>%     #alternate code: sir_function(num, denom)
  unnest_wider(result)# %>% 
# mutate(across(c(midp, sir, llimit, ulimit), ~ round(., 3)))    #format to 3 decimal places

head(df1)




