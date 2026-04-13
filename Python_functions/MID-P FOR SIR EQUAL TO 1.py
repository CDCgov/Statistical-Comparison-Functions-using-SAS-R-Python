# May 2025
# Compare SIR from 1 or other nominal values
# SAS Macro downloaded from https://www.cdc.gov/nhsn/ps-analysis-resources/index.html
# Output: p-value & 95% CI
# Converting from SAS Macro "SIRcomp" to Python function by Qunna Li
# Tested by Rob Wild


import math
import scipy.stats as stats
import numpy as np
import pandas as pd

# The function was developed using the following versions
# scipy version: 1.4.1
# pandas version: 1.0.1
# numpy version: 1.18.1

def SIRcomp(OBS, EXP):
    """
    Usage: Two-sided p-value for SIR compared to 1 or other nominal value.
    
    Method:
        Exact p-value with mid-p method based on the Poisson distribution will be computed only
        when observed numbers of events are <=100 (point estimate of mu <=100).
        Otherwise, Byar approximation method (based on the large sample method) is used for computational efficiency,
        assuming similar results with mid-P method. 
        The one-tailed P-values thus obtained are doubled to provide two-tailed P-values in the final results.

    How to use the function?
        A facility’s SIR can be compared to a single nominal value, such as a “goal” or “target” SIR, using this function. 
        To do this, you will need to re-calibrate the number of predicted events for each hospital based on the target SIR. 

        For example, facility A has '2' observed events and '2.04' predicted events, 
        resulting in the SIR of 0.98 for a given time period,
        and we would like to compare this SIR value (0.98) to the target SIR of 0.80. 
        - First, we must multiply the number of predicted events (2.04) by the target SIR (0.80), 
         and obtain the newly calculated number of predicted events (2.04 x 0.80 = 1.63).
        - then, for parameter values in the function, the number of observed events is '2' for 'OBS' 
         and number of predicted events is '1.63' for 'EXP'.
        - Last, apply the function to obtain the p-value. If the two-sided p-value is <= 0.05, 
         we can conclude that facility A’s SIR of 0.98 is significantly different from the target SIR of 0.80.
  
    Reference: 
        1. JH ABRAMSON. WINPEPI PROGRAMS. DESCRIBE MANUAL (VERSION 2.42), PAGE-52. Available at 'http://www.brixtonhealth.com/pepi4windows.html'
        2. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com.
        3. Rothman KJ, Boice JD Jr: Epidemiologic analysis with a programmable calculator. NIH Pub No. 79-1649. Bethesda, MD: National Institutes of Health, 1979;31-32.
        4. Rothman KJ, Greenland S.  Modern Epidemiology, 2nd Edition.  Lippincott-Raven Publishers, Philadelphia, 1998.
        5. GEOFFREY RC, SHU-YING Y. MID-P CONFIDENCE INTERVALS FOR THE POISSON EXPECTATION. STATISTICS IN MEDICINE, 13,2189-2203 (1994)

    Output: the function returns SIR_pval, SIR, SIR_L, SIR_U
        Parameters:
            OBS: number of observed events
            EXP: number of predicted events

    Note: No results will be generated if entry values are missing or impossible values.
    Note: When SIR==0, the lower limit of 95%CI (SIR_L) is set to missing by default. A user may change it to SIR_L==0 if desired.
    """
    e = math.e
    
    # Ensure result always contains the expected keys
    result = {"SIR_pval": np.nan, "SIR": np.nan, "SIR_L": np.nan, "SIR_U": np.nan}

    # Handling missing data and cases where EXP < 1
    if pd.isna(OBS) or OBS<0 or pd.isna(EXP) or EXP < 1 or (OBS % 1 !=0):
        return result  # Ensuring a uniform structure

    # Case when OBS <= 100
    elif OBS <= 100:
        total = 0
        OBS=int(OBS)
        if OBS < EXP:
            total = sum((e**-EXP) * (EXP**k) / math.factorial(k) for k in range(OBS))
            num = (e**-EXP) * (EXP**OBS)
            aa = (num / math.factorial(OBS)) * 0.5
            result["SIR_pval"] = 2 * (aa + total)
        else:
            total = sum((e**-EXP) * (EXP**k) / math.factorial(k) for k in range(OBS))
            num = (e**-EXP) * (EXP**OBS)
            aa = (num / math.factorial(OBS)) * 0.5
            result["SIR_pval"] = 2 * (1 - (aa + total))

        # Ensure valid range
        result["SIR_pval"] = min(max(result["SIR_pval"], 0), 1)

    # Case when OBS > 100, using Byar Poisson approximation
    else:
        OBS=int(OBS)
        OBS_ = OBS if OBS > EXP else OBS + 1
        prob_norm = stats.norm.cdf(math.sqrt(9 * OBS_) * (1 - 1 / (9 * OBS_) - (EXP / OBS_) ** (1 / 3)))

        result["SIR_pval"] = 2 * prob_norm if OBS <= EXP else 2 * (1 - prob_norm)
        result["SIR_pval"] = min(max(result["SIR_pval"], 0), 1)

    # 95% CI Calculation
    def fish(z, x1, x2):
        q, tot, s, k = 1, 0, 0, 0
        while k < z or q > (tot * 10**-10):
            tot += q
            if x1 <= k <= x2:
                s += q
            if tot > 10**30:
                s /= 10**30
                tot /= 10**30
                q /= 10**30
            k += 1
            q *= z / k
        return s / tot

    def fish2(z, x1, x2):
        return fish(z, x1, x2) * 0.5

    if EXP >= 1:
        result["SIR"] = OBS / EXP

        # Lower tail
        v, dv, p = 0.5, 0.5, 2.5 / 100
        while dv > 10**-5:
            dv /= 2
            poisp = fish((1 + OBS) * v / (1 - v), OBS + 1, 10**10)
            poisp2 = fish2((1 + OBS) * v / (1 - v), OBS, OBS)
            v = v - dv if (poisp + poisp2) > p else v + dv

        result["SIR_L"] = ((1 + OBS) * v / (1 - v)) / EXP if OBS != 0 else np.nan

        # Upper tail
        v, dv, p = 0.5, 0.5, 2.5 / 100
        while dv > 10**-5:
            dv /= 2
            poisp = fish((1 + OBS) * v / (1 - v), 0, OBS - 1)
            poisp2 = fish2((1 + OBS) * v / (1 - v), OBS, OBS)
            v = v - dv if (poisp + poisp2) < p else v + dv

        result["SIR_U"] = ((1 + OBS) * v / (1 - v)) / EXP

    return result


##################################################################################
# The following code used for testing the function


# Create a sample dataframe

df = pd.DataFrame({'obs': [2, np.nan, 3, 1,0,1, 105,110, -1, 1], 'pred': [4, 5, np.nan, 0.8,1,0, 101, 115,1, -2]})

# Apply SIRcomp function to each row
df[['SIR_pval','SIR', 'SIR_L', 'SIR_U']] = df.apply(lambda row: pd.Series(SIRcomp(row['obs'], row['pred'])), axis=1)

# Round values if needed
columns_to_round_list=['SIR_pval','SIR', 'SIR_L', 'SIR_U']

for i in columns_to_round_list:
    if i == 'SIR_pval':
        df[i]=df[i].round(4)
    else:
        df[i]=df[i].round(3)
  
print(df)


