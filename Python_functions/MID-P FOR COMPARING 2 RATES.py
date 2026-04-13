# May 2025
# Converting SAS Macro "IDRcomp" to Python code
# Compare two Incidence Density Rates
# SAS Macro downloaded from https://www.cdc.gov/nhsn/ps-analysis-resources/index.html
# Convert from sas macro "IDRcomp" to Python function "IDRcomp" by Qunna Li
# Tested by Rob Wild



import scipy
import scipy.stats as stats
import numpy as np
import pandas as pd

# The function was developed using the following versions
# scipy version: 1.4.1
# pandas version: 1.0.1
# numpy version: 1.18.1

def two_rates(o1, pt1, o2, pt2):
    """
    Usage: compare two incidence density rates. 
    Method: this method is based on the fact that distribution of two poisson variates conditional on their sum is binomial.
    Reference: 
        1. Erich L. Lehmann, Joseph P. Romano. Testing Statistical Hypotheses. Wiley, New York, 1959(Third edition, 2005, Springer Texts in Statistics) 
        2. John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
        3. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com, updated 2013/03/21
        4. Austin, Harland (Emory). Epidemiology method-II: Statistical issues for Density type follow-up studies.

    Input:         
            o1: observed events in group-1
            pt1: person-time in group-1
            o2: observed events in group-2
            pt2: person-time in group-2

    Output: 
            MID_P: significance test for comparing between 2 rates
            RATE_RATIO: Rate in group-2 / Rate in group-1
            LL: lower limit of 95%CI of rate ratio
            UL: Upper limit of 95%CI of rate ratio
    Note:   When RATE_RATIO=0, the lower limit (LL) of 95%CI is set to missing by default. A user may change it to LL=0 if desired.   

    """
    def bin_p(x1, x2, pp, n):
        """Statistical function for interval calculation"""
        q = pp / (1 - pp)
        k, v, s, tot = 0, 1, 0, 0

        while k <= n:
            tot += v
            if x1 <= k <= x2:
                s += v
            if tot > 10**30:
                s /= 10**30
                tot /= 10**30
                v /= 10**30
            k += 1
            v = v * q * (n + 1 - k) / k

        return s / tot

    def bin_p2(x1, x2, pp, n):
        """Statistical function for interval calculation"""
        return bin_p(x1, x2, pp, n) * 0.5
     # Ensure result always contains the expected keys
    result = {"MID_P": np.nan, "RATE_RATIO": np.nan, "LL": np.nan, "UL": np.nan}   
        
    if pd.isna(o1) or o1<0 or pd.isna(pt1) or pt1 <= 0 or pd.isna(o2) or o2<0 or pd.isna(pt2) or pt2 <= 0 or (o1==0 and o2==0) or o1>pt1 or o2>pt2 or (o1 % 1 !=0) or (o2 % 1 !=0) or (pt1 % 1 !=0) or (pt2 % 1 !=0):
        return result   # Ensuring a uniform structure

    else:
        vn = o1 + o2
        vp = o2 / vn

        # Calculate Mid-P 
        ratio1 = o1 / pt1
        ratio2 = o2 / pt2
        total = pt1 + pt2
        o = o1 + o2

        # Take the larger of RISKs (Testing for Ha: RR>1)
        if ratio1 >= ratio2:
            p1 = 1 - stats.binom.cdf(o1, o, pt1 / total)
            p2 = 0.5 * (stats.binom.cdf(o1, o, pt1 / total) - stats.binom.cdf(o1 - 1, o, pt1 / total))
        else:
            p1 = 1 - stats.binom.cdf(o2, o, pt2 / total)
            p2 = 0.5 * (stats.binom.cdf(o2, o, pt2 / total) - stats.binom.cdf(o2 - 1, o, pt2 / total))

        mid_p = 2 * (p1 + p2)
        # Ensure MID-P within 0 to 1
        mid_p = min(max(mid_p, 0), 1)

        result["MID_P"] = mid_p 

         # Rate Ratio Calculation
        if o1 == 0:
            result.update({"RATE_RATIO": np.nan, "LL": np.nan, "UL": np.nan})
        else:
            rate_ratio = (o2 / pt2) / (o1 / pt1)

          # Confidence Interval Calculation
            tol = 1e-5

            # Lower Bound (DL)
            if o2 == 0:
                dl = 0
            else:
                p2 = vp / 2
                low, high = 0, vp
                while (high - low) > tol:
                    if bin_p(o2 + 1, vn, p2, vn) + bin_p2(o2, o2, p2, vn) > 0.025:
                        high = p2
                        p2=(low+p2)/2
                    else:
                        low = p2
                        p2=(high+p2)/2
                dl = p2

            # Upper Bound (DU)
            if o2 == vn:
                du = 1
            else:
                p3 = (1 + vp) / 2
                low, high = vp, 1
                while (high - low) > tol:
                    if bin_p(0, o2 - 1, p3, vn) + bin_p2(o2, o2, p3, vn) < 0.025:
                        high = p3
                        p3=(low+p3)/2
                    else:
                        low = p3
                        p3=(high+p3)/2
                du = p3

            ll = (dl * pt1) / ((1 - dl) * pt2) if (not pd.isna(dl) and dl !=0) else np.nan
            ul = (du * pt1) / ((1 - du) * pt2) if not pd.isna(du) else np.nan

            result.update({"RATE_RATIO": rate_ratio, "LL": ll, "UL": ul})
   
    
    return result


# Following code used for testing the function

# Sample DataFrame
data = pd.DataFrame({
    "O1": [15, 20, np.nan,30,  4, 0, 1,1, 1, 1, 5],  # Observed events in group 1
    "PT1": [1000, 1200, 1500, 900,39, 10, 50,np.nan,0, 10, 4],  # Person-time for group 1
    "O2": [25, 30, 45, 10,np.nan, 1, 0, 1, 1, 1, 5],  # Observed events in group 2
    "PT2": [1200, 1400, 1800, np.nan,70, 20, 30, 20, 10, 0, 5]  # Person-time for group 2
})


# Apply the function to each row
data[['MID_P', 'RATE_RATIO', 'LL', 'UL']] = data.apply(lambda row:pd.Series(two_rates(row["O1"], row["PT1"], row["O2"], row["PT2"])), axis=1)

print(data)


