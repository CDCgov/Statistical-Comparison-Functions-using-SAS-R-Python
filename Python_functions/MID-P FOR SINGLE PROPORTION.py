# June 2025
# Function: calculate 2 SIDED P-VALUE FOR SINGLE PROPORTION
# SAS Macro downloaded from https://www.cdc.gov/nhsn/ps-analysis-resources/index.html
# Convert from sas macro to Python function by Qunna Li
# Tested by Rob Wild

import numpy as np
import pandas as pd



# The function was developed using the following versions
# pandas version: 1.0.1
# numpy version: 1.18.1

def binom_ci(vX, vN):
    """
    Usage: calculate 2 SIDED P-VALUE FOR a SINGLE PROPORTION. 
    Reference: 
        1. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version 2.3.1. www.OpenEpi.com
        2. Rothman KJ,  Boice JD  Jr:  Epidemiologic analysis with a programmable calculator. NIH Pub No. 79-1649.  Bethesda, MD:  National Institutes of Health, 1979;31-32.
        3. Fleiss JL.  Statistical Methods for Rates and Proportions, 2nd Ed. John Wiley & Sons, New York, 1981.
        4. Newcombe RG.  Two-sided confidence intervals for the single proportion: comparison of seven methods.  Statistics in Medicine 1988;17:857-872.
        5. Stein Vollset. CONFIDENCE INTERVALS FOR A BINOMIAL PROPORTION STATISTICS IN MEDICINE,12,809-824 (1993)*/

    Input parameters:
        vX = numerator
        vN = denominator
        
    Output: 
        vp = proportion
        dl = 95% lower bound of proportion
        du = 95% upper bound of proportion

    Note: When vp==0, the 95% lower bound of proportion (dl) is set to missing by default. A user may change it to dl==0 if desired.
     """
    
    def BinP(x1, x2):
        q = p2 / (1 - p2)
        k, v, s, tot = 0, 1, 0, 0
        while k <= vN:
            tot += v
            if x1 <= k <= x2:
                s += v
            if tot > 10**30:
                s /= 10**30
                tot /= 10**30
                v /= 10**30
            k += 1
            v = v * q * (vN + 1 - k) / k
        return s / tot

    def BinP2(x1, x2):
        return BinP(x1, x2) * 0.5

    # Handling missing values and illogical values
    result = {"vp": np.nan, "dl": np.nan, "du": np.nan}
    if pd.isna(vX) or pd.isna(vN) or vN==0 or vX < 0 or vN < 0 or vN < vX or (vX % 1 != 0) or (vN % 1 != 0):
        return result
    else:
        vp = vX / vN
 
        # Lower bound calculation
        if vX > 0:
            p2 = vp / 2
            vsL, vsH, p = 0, vp, 2.5 / 100
            while (vsH - vsL) > 1e-5:
                if (BinP(vX+1, vN) + BinP2(vX, vX)) > p:
                    vsH = p2
                    p2 = (vsL + p2) / 2
                else:
                    vsL = p2
                    p2 = (vsH + p2) / 2
            dl = p2
        else:
            dl=np.nan

        # Upper bound calculation
        if vX == vN:
            du=1
        else:
            p2 = (1 + vp) / 2
            vsL, vsH, p = vp, 1, 2.5 / 100
            while (vsH - vsL) > 1e-5:
                if (BinP(0, vX - 1) + BinP2(vX, vX)) < p:
                    vsH = p2
                    p2 = (vsL + p2) / 2
                else:
                    vsL = p2
                    p2 = (vsH + p2) / 2
            du = p2

        result.update({"vp": vp, "dl": dl, "du": du})
    return result

##########################################################################################################
# Apply the function to a test dataframe

# Sample DataFrame
data = pd.DataFrame({
    "num": [5, 0, np.nan, -1, 2.3, 5, 5,5, 5, 0],  # Numerator
    "denom": [10, 5, 5, 5,5, 4, np.nan, 5.5, -10, 0]  # Denominator
})


# Apply the function to each row
data[['vp', 'dl', 'du']] = data.apply(lambda row:pd.Series(binom_ci(row["num"], row["denom"])), axis=1)

print(data)