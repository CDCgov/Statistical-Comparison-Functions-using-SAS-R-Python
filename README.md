# Statistical-Comparison-Functions-using-SAS-R-Python
Organization: CDC/NCEZID/DHQP/SB

Contact email: ncezid_shareit@cdc.gov

Description: The Statistics Team members within the Surveillance Branch (SB), DHQP, have developed a series of statistical comparison functions. These functions were developed to conduct statistical tests for epidemiological measures and summary statistics, especially those consistent with the methods used by National Healthcare Safety Network (NHSN).

Initial statistical comparison functions test several common epidemiology measures such as incidence density rates, proportions, and ratios such as standardized metrics known as the standardized infection ratio (SIR), standardized utilization ratio (SUR), standardized antimicrobial administration ratio (SAAR), etc., and all statistical tests are based on the mid-p exact method. These functions were originally written using the SAS macro language in early 2010s. Six SAS macros were distributed for internal analysis use within CDC’s SB, and were integrated in the NHSN web application to support users, as well as published on the NHSN website for public use (https://www.cdc.gov/nhsn/ps-analysis-resources/index.html). In early 2026, NHSN updated functionality of these six macros and developed equivalent Python and R versions comparable to the original SAS implementations.

To further complement the statistical analysis of epidemiological data, an additional function written in R was developed to compare summary statistics such as the medians and empirical cumulative distribution functions (ECDFs) among two groups for any numeric variable. This function uses nonparametric methods such as Brown-Mood median test and both the Anderson-Darling and Dowd DTS tests for comparing ECDFs.

Subfolders for SAS, Python and R have been created to organize the statistical functions by corresponding programming languages and quick reference guide containing example use cases are also provided in respective subfolders.

We’ve chosen to use the term ‘function’ to more generically refer and describe these routines although they were originally developed as SAS macros at first.
Languages: R, Python, Other
