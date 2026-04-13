# June 2025
# Function: COMPARE 2 PROPORTIONS (2 X 2 TABLE) by mid-P exact method
# SAS Macro downloaded from https://www.cdc.gov/nhsn/ps-analysis-resources/index.html
# Convert from sas macro to Python function by Qunna Li
# Tested by Rob Wild


import math
import numpy as np
import pandas as pd

# The function was developed using the following versions
# pandas version: 1.0.1
# numpy version: 1.18.1

def TWOBYTWO(A, B, C, D):
    """
    Usage: Compare 2 proportions (2 X 2 table) by mid-P exact method.
    Method: Hypothesis testing of comparing 2 proportions by Mid-p method is based on hypergeometric distribution.
    Use case scenarios: 
      Analyzing 2 x 2 contingency table that shows the findings in two independent groups(within a single stratum).
      The data may be derived from an observational study (cross-sectional or cohort) or from a trial.
    Reference:
      1.David G. Kleinbaum, Lawrence L. Kupper. Epidemiologic Research: Principles and Quantitative Methods.
      2.Dean AG, Sullivan KM, Soe Minn Minn. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. 
      www.OpenEpi.com, updated 2013/03/21
      3.John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
      4.Bernard Rosner. Fundamentals of Biostatistics' (5th edition) (Example 10.20, page 375). 
      Two-tailed p-value calculated as 2 times whichever is smallest: left-tail, right-tail, or 0.5. 
      It tends to agree closely with Yates Chi-Square p-value.
      5.David Martin, Harland Austin. An Efficient Program for Computing Conditional Maximum Likelihood Estimates and
      Exact Confidence Limits for a Common Odds Ratio. Epidemiology Vol. 2, No. 5 (Sep., 1991), pp. 359-362
      
    Input parameters:
      A= E+,D+
      B= E-,D+
      C= E+,D-
      D= E-,D-
      
      
                                          Disease/Outcome 
                      --------------------------------------------------------
                                 |    Yes                No           |                
                      -----------|------------------------------------|-------         
               Exposure     Yes  |     A                  C           |                
                            No   |     B                  D           |                
                      -----------|------------------------------------|-------         
    
    Output: midP = TWOTAILED MID-P value 
    """
    
    def LNfact(z):
        f, lnfact = 0, 0
        if z < 2:
            lnfact = 0
        elif z < 17:
            f = z
            while z > 2:
                z -= 1
                f *= z
            lnfact = math.log(f)
        else:
            ln_pi2 = math.log(2 * math.pi)
            lnfact = (z + 0.5) * math.log(z) - z + ln_pi2 / 2 + 1 / (12 * z) - 1 / (360 * (z**3)) + 1 / (1260 * (z**5)) - 1 / (1680 * (z**7))

        return lnfact

    result = {"midP": np.nan}
    
    if (A == 0 and B == 0 and C == 0 and D == 0) or pd.isna(A) or pd.isna(B) or pd.isna(C) or pd.isna(D) or  \
              A < 0 or B < 0 or C < 0 or D < 0 or (A % 1 != 0) or (B % 1 != 0) or (C % 1 != 0) or (D % 1 != 0):
        return result
    
    else:
        Cell_A, Cell_B, Cell_C, Cell_D = A, B, C, D
        Cell_r1 = Cell_A + Cell_B
        Cell_r2 = Cell_C + Cell_D
        Cell_c1 = Cell_A + Cell_C
        Cell_c2 = Cell_B + Cell_D
        t = Cell_A + Cell_B + Cell_C + Cell_D

        lo_slop = min(Cell_A, Cell_D)
        hi_slop = min(Cell_B, Cell_C)

        ln_prob1 = LNfact(Cell_r1) + LNfact(Cell_r2) + LNfact(Cell_c1) + LNfact(Cell_c2) - LNfact(t)

        fisher_p, left_p, left_p1, right_p, right_p1 = 0, 0, 0, 0, 0
        k = Cell_A - lo_slop

        while k <= Cell_A + hi_slop:
            p = math.exp(ln_prob1 - LNfact(k) - LNfact(Cell_r1 - k) - LNfact(Cell_c1 - k) - LNfact(k + Cell_r2 - Cell_c1))

            if k <= Cell_A:
                left_p += p
            if k < Cell_A:
                left_p1 += p
            midp_left = ((left_p - left_p1) * 0.5 + left_p1)

            if k >= Cell_A:
                right_p += p
            if k > Cell_A:
                right_p1 += p
            midp_right = ((right_p - right_p1) * 0.5 + right_p1)

            k += 1

        fisher_p = 2 * min(left_p, right_p)
        mid_p = min(max(2 * min(midp_left, midp_right), 0), 1)
        result.update({"midP": mid_p})

    return result


##########################################################################################################
# Apply the function to a test dataframe


# sample dataframe
data=pd.DataFrame({
    "A": [-1, 10, 1, 1, 1.6, 0, np.nan, 10, 10, 10, 5, 6, 12, 10, 12, 2, 3, 25, 100, 0],
    "B": [20, 0, np.nan, 2, 4, 5, 5, -2, 1.6, 5, 5, 8, 20, 23, 22, 5, 6, 50, 500, 0],
    "C": [11, 7, 3, 3, 3, 7, 8, 3, 3, np.nan, -3, 3.6, 0, 2, 23, 6, 7, 100, 220, 0],
    "D": [21, 20, 4, 4, 6, 9, 9, 8, 9, 5, 5, 6, 6, np.nan, 20.5, -1, 0, 150, 600, 0]
})


# Apply the function to each row
data['midP']=data.apply(lambda row:pd.Series(TWOBYTWO(row["A"], row["B"], row["C"], row["D"])), axis=1)

print(data)

