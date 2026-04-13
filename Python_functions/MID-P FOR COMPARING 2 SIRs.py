
# May 2025
# Compare 2 standardized ratio (2 SIRs) based on mid-p exact method
# SAS Macro downloaded from https://www.cdc.gov/nhsn/ps-analysis-resources/index.html
# Convert from sas macro to Python function by Qunna Li
# Tested by Rob Wild


import scipy
import scipy.stats as stats
import numpy as np
import pandas as pd

# The function was developed using the following versions
# scipy version: 1.4.1
# pandas version: 1.0.1
# numpy version: 1.18.1

def binom(o1, e1, o2, e2):
    """
    Usage: compare two standardized ratio (two SIRs). 
    Method: based on the assumption that the distribution of two poisson variates conditional on their sum is binomial.
    Reference: 
        1. Breslow NE, Day NE. Statistical methods in cancer research. Volume II--The design and analysis of cohort studies. IARC Sci Publ. 1987;(82):1-406. PMID: 3329634.
        http://www.iarc.fr/en/publications/pdfs-online/stat/sp82/SP82_vol2-3.pdf
        2. Erich L. Lehmann, Joseph P. Romano. Testing Statistical Hypotheses. Wiley, New York, 1959 (Third edition, 2005, Springer Texts in Statistics) 
        3. John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
        4. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com, updated 2013/03/21

    Input parameters:
        o1= observed events in group1
        e1= predicted events in group1
        o2= observed events in group2
        e2= predicted events in group2
        
    Output: 
        midP= TWOTAILED MID-P value
        RATIO= SIR in group2 / SIR in group1 = SIR2/SIR1 = (o2/e2)/(o1/e1)
        LL= 95% lower bound of RATIO
        UL= 95% Upper bound of RATIO 

    Note: When RATIO=0, the lower limit (LL) of 95%CI is set to missing by default. A user may change it to LL=0 if desired.      
    """
    
            
    def binp(x1, x2, pp, n):
        q = pp / (1 - pp)
        k,v,s,tot = 0, 1, 0, 0

        while k <= n:
            tot += v
            if x1 <= k <= x2:
                s += v
            if tot > 10**30:
                s /= 10**30
                tot /= 10**30
                v /= 10**30
            k += 1
            v *= q * (n + 1 - k) / k

        return s / tot

    def binp2(x1, x2, pp, n):
        return binp(x1, x2, pp, n) * 0.5


    result = {"midP": np.nan, "RATIO": np.nan, "LL": np.nan, "UL": np.nan}
    
    if (o1 == 0 and o2 == 0) or pd.isna(o1) or pd.isna(e1) or pd.isna(o2) or pd.isna(e2) or e1<1 or e2<1 or o1<0 or o2<0 or (o1 % 1 !=0) or (o2 % 1 !=0):
        return result

    else:
        vn = o1 + o2
        vp = o2 / vn

        # midP value for hypothesis testing
        ratio1 = o1 / e1
        ratio2 = o2 / e2
        o = o1 + o2
        t = e1 + e2

        if ratio1 >= ratio2:
            p1 = 1 - stats.binom.cdf(o1, o, e1/t)
            p2 = 0.5 * (stats.binom.cdf(o1, o, e1/t) - stats.binom.cdf(o1-1, o, e1/t))
        else:
            p1 = 1 - stats.binom.cdf(o2, o, e2/t)
            p2 = 0.5 * (stats.binom.cdf(o2, o, e2/t) - stats.binom.cdf(o2-1, o, e2/t))

        result['midP']=2 * (p1 + p2)
        
        # Ensure valid range
        result["midP"] = min(max(result['midP'],0), 1)

        
        # Ratio and interval estimation (SIR2/SIR1, apriori assumption: SIR2>SIR1)

        if o1 == 0:
            result.update({"RATIO": np.nan, "LL": np.nan, "UL": np.nan})
        else:
            ratio = (o2/e2) / (o1/e1)

            if o2 == 0:
                dl = 0
            else:
                p2 = vp / 2
                vsL = 0
                vsH = vp
                p = 2.5 / 100
                while (vsH - vsL) > 10**-5:
                    if binp(o2+1, vn, p2, vn) + binp2(o2, o2, p2, vn) > p:
                        vsH = p2
                        p2 = (vsL + p2) / 2
                    else:
                        vsL = p2
                        p2 = (vsH + p2) / 2
                dl = p2

            if o2 == vn:
                du = 1
            else:
                p3 = (1 + vp) / 2
                vsL = vp
                vsH = 1
                p = 2.5 / 100
                while (vsH - vsL) > 10**-5:
                    if binp(0, o2-1, p3, vn) + binp2(o2, o2, p3, vn) < p:
                        vsH = p3
                        p3 = (vsL + p3) / 2
                    else:
                        vsL = p3
                        p3 = (vsH + p3) / 2
                du = p3

            ll = (dl * e1) / ((1 - dl) * e2)
            ul = (du * e1) / ((1 - du) * e2)

            if ratio == 0:
                ll = np.nan

            result.update({"RATIO": ratio, "LL": ll, "UL": ul})

    return result



##########################################################################################################
# Apply the function to a test dataframe


# Sample DataFrame
data = pd.DataFrame({
    "o1": [1.1, 1.0, 3.0, 10.0, 30, 2, 2, -2,3,3,np.nan,2,2,0,1,0,1],  # Observed events in group 1
    "e1": [10,4,2,100,20, 4,4,4,4,0.6,4,np.nan, 4,4,2,5,2],  # predicted events in group1
    "o2": [4,2.1,4,40,40,3,-3,3,3,3,3,3,3,1,0,0,np.nan],  # Observed events in group 2
    "e2": [5,5,4,50,40,-5,5,5,0.5,5,5,5,np.nan,5,5,9,20]  # predicted events in group2
})


# Apply the function to each row
data[['midP', 'RATIO', 'LL', 'UL']] = data.apply(lambda row:pd.Series(binom(row["o1"], row["e1"], row["o2"], row["e2"])), axis=1)

print(data)

