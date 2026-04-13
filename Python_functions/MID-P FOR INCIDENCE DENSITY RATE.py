
# June 2025
# Calculate 95%CI of Incidence Density Rate by mid-P exact method
# SAS Macro downloaded from https://www.cdc.gov/nhsn/ps-analysis-resources/index.html
# Convert from sas macro to Python function by Qunna Li
# Tested by Rob Wild


import numpy as np
import pandas as pd

# The function was developed using the following versions
# pandas version: 1.0.1
# numpy version: 1.18.1

def rateCIComp(numer, denom):
    """
    Usage: Calculate Incidence Density Rate and its 95%CI
    Method: mid-P exact method
    Reference:
    1. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version 2.3.1. 
    www.OpenEpi.com
    2. Rothman KJ,  Boice JD  Jr:  Epidemiologic analysis with a programmable calculator.
    NIH Pub No. 79-1649.  Bethesda, MD:  National Institutes of Health, 1979;31-32.
    3. Rothman KJ, Greenland S.  Modern Epidemiology, 2nd Edition.  Lippincott-Raven Publishers, Philadelphia, 1998.
    
    Input parameters:
        numer: numerator
        denom: denominator
        
    Output:
        rate: incidence density rate per 1000
        rate_l: lower limit of 95%CI
        rate_u: upper limit of 95%CI
        
    Note: When rate=0, the lower limit (rate_l) of 95%CI is set to missing by default. 
        A user may change it to rate_l=0 if desired.
        
    """

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

    # Handling of missing/impossible values and relation
    result = {"rate": np.nan, "rate_l": np.nan, "rate_u": np.nan}
    if (
        (numer == 0 and denom == 0) or
        pd.isna(numer) or pd.isna(denom) or
        numer < 0 or denom <= 0 or
        numer > denom or
        (numer % 1 != 0) or (denom % 1 != 0)
    ):
        return result

    else:
        # Crude rate calculation
        rate = (numer / denom) * 1000

        # Lower CI bound
        if numer == 0:
            rate_l=np.nan
        else:
            v, dv, p = 0.5, 0.5, 2.5 / 100
            while dv > 10**-5:
                dv /= 2
                z = (1 + numer) * v / (1 - v)
                poisP = fish(z, numer + 1, 10**10)
                poisP2 = fish2(z, numer, numer)
                if (poisP + poisP2) > p:
                    v -= dv
                else:
                    v += dv
            rate_l = ((1 + numer) * v / (1 - v)) / denom * 1000

        # Upper CI bound
        v, dv, p = 0.5, 0.5, 2.5 / 100
        while dv > 10**-5:
            dv /= 2
            z = (1 + numer) * v / (1 - v)
            poisP = fish(z, 0, numer - 1)
            poisP2 = fish2(z, numer, numer)
            if (poisP + poisP2) < p:
                v -= dv
            else:
                v += dv
        rate_u = ((1 + numer) * v / (1 - v)) / denom * 1000
        
        result.update({"rate": rate, "rate_l": rate_l, "rate_u": rate_u})

    return result


##########################################################################################################
# Apply the function to a test dataframe

# Sample DataFrame
data = pd.DataFrame({
    "numer": [0, 10, 1.1, 1, np.nan, 0, 0, -4, 4, 4, 20, 20],
    "denom": [10, 10, 10, 4, 2, 0, 10, 50, np.nan, -50, 15, 50.5]
})


# Apply the function to each row
data[['rate', 'rate_l', 'rate_u']] = data.apply(lambda row:pd.Series(rateCIComp(row["numer"], row["denom"])), axis=1)

print(data)

