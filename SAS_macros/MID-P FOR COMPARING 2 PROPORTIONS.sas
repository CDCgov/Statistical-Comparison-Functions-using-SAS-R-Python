/*********************************************************************************************************************

						COMPARE 2 PROPORTIONS (2 X 2 TABLE) by mid-P exact method

*********************************************************************************************************************
Function: This macro can be used to analyse any simple 2 x 2 contingency table that shows the findings in two independent groups(within a single stratum). 
The data may be derived from an observational study (cross-sectional or cohort) or from a trial. 

Use case scenarios:
For hypothesis testing of comparing 2 proportions by Mid-p method based on hypergeometric distribution.

Developed by Minn M. Soe(Lead Epidemiologist, Surveillance Branch, NHSN, CDC), May 2025

Validation: with OpenEpi, WinPEPI, www.Statpages.com

**********************************************************************************************************************
Data input layout:
A= E+,D+
B= E-,D+
C= E+,D-
D= E-,D-
   								 Disease/Outcome								   
   			             |    Yes                No           |                
			  -----------|------------------------------------|-------         
       Exposure 	Yes	 |     A                  C			  |                
  	   				No   |     B                  D			  |                
  			  -----------|------------------------------------|-------         

Output:
-midp(2 sided)

**********************************************************************************************************************

Reference: 
1.David G. Kleinbaum, Lawrence L. Kupper. Epidemiologic Research: Principles and Quantitative Methods. 
2.Dean AG, Sullivan KM, Soe Minn Minn. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com, updated 2013/03/21
3.John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
4.Bernard Rosner. Fundamentals of Biostatistics' (5th edition) (Example 10.20, page 375). Two-tailed p-value calculated 
as 2 times whichever is smallest: left-tail, right-tail, or 0.5. It tends to agree closely with Yates Chi-Square p-value.
5.David Martin, Harland Austin. An Efficient Program for Computing Conditional Maximum Likelihood Estimates and 
Exact Confidence Limits for a Common Odds Ratio. Epidemiology Vol. 2, No. 5 (Sep., 1991), pp. 359-362 
**********************************************************************************************************************/;


%MACRO TWOBYTWO(A=,B=,C=,D=);

%MACRO LNfact(Z=);
F=0;LNFACT=0;z=&z;
if(z<2) THEN LNFACT=0;

ELSE if(z<17) THEN DO;
 	f=z; 
	DO while(z>2);  z=z-1; f=f*z ; 
	END;
	LNFACT=LOG(f) ;
END;

ELSE DO;
  LNFACT=(z+0.5)*LOG(z) - z + LnPi2/2 + 1/(12*z) - 1/(360*(z**3)) + 1/(1260*(z**5)) - 1/(1680*(z**7))   ;
END;
%MEND;



***********INITIAL SETTING;
*HANDLING OF MISSING VALUES;
IF  (&A = 0 AND &B = 0 AND &C =0 AND &D = 0) or &A=. or &B=. or &C=. or &D=. or &A<0 or &B<0 or &C<0 or &D<0  or mod(&A, 1) ^= 0  or mod(&B, 1) ^= 0  or mod(&C, 1) ^= 0  or mod(&D, 1) ^= 0  THEN DO;
MID_P=.;
END;

ELSE DO;

Cell_A=0; Cell_B=0; Cell_C=0;  Cell_D=0;
Cell_r1=0;   Cell_r2=0;   Cell_c1=0;   Cell_c2=0;   t=0;


*INPUT PARAMETERS;
  Cell_A = &A;
  Cell_B = &B;
  Cell_C = &C;
  Cell_D = &D;
  Cell_r1 = Cell_A+Cell_B;
  Cell_r2 = Cell_C+Cell_D;
  Cell_c1 = Cell_A+Cell_C;
  Cell_c2 = Cell_B+Cell_D;
  t = Cell_A+Cell_B+Cell_C+Cell_D;

*SETTING THE DEFAULT VALUES;
Pi=3.141592653589793;  Pi2=2*Pi; LnPi2 = LOG(Pi2);
E1=0; E2=0; E3=0; E4=0; E5=0; 


***************************HYPOTHESIS TESTING;
*STAT FUNCTIONS;
*function CalcStats(form);
  LoSlop = Cell_A; if(Cell_D<Cell_A) THEN  LoSlop = Cell_D;
  HiSlop = Cell_B; if(Cell_C<Cell_B) THEN  HiSlop = Cell_C;
  *LnProb1 = %LnFact(Cell_r1) + %LnFact(Cell_r2) + %LnFact(Cell_c1) + %LnFact(Cell_c2) - %LnFact(t);
  %LnFact(Z=Cell_r1);E1=LNFACT;  %LnFact(Z=Cell_r2);E2=LNFACT;  %LnFact(Z=Cell_c1);E3=LNFACT; 
%LnFact(Z=Cell_c2);E4=LNFACT;  %LnFact(Z=t);E5=LNFACT; 
  LnProb1 = E1+E2+E3+E4 -E5;

  *SingleP = Exp( LnProb1 - %LnFact(Cell_A) - %LnFact(Cell_B) - %LnFact(Cell_C) - %LnFact(Cell_D) );
  FisherP=0;  LeftP=0;  LeftP1=0; RightP=0;  RightP1=0; SumCheck=0;  
  k = Cell_A - LoSlop;

 *EXACT TESTS;
  DO while( k<=Cell_A+HiSlop );
    %LnFact(Z=k); E1=LNFACT;  %LnFact(Z=Cell_r1-k) ; E2=LNFACT;   %LnFact(Z=Cell_c1-k); E3=LNFACT; 
	%LnFact(Z=k+Cell_r2-Cell_c1);E4=LNFACT; 
	P = Exp( LnProb1 -E1 -E2 -E3 -E4);

    *SumCheck = SumCheck + P;
    if( k<=Cell_A  ) THEN LeftP   = LeftP   + P ;
	if( k<Cell_A  ) THEN LeftP1   = LeftP1   + P ;
	midp_left=((LEFTP - LEFTP1)*0.5 + LEFTP1);
	
	if( k>=Cell_A  ) THEN RightP  = RightP  + P ;
    if( k>Cell_A  ) THEN RightP1   = RightP1   + P ;
	midp_Right=((RightP - RightP1)*0.5 + RightP1);
	
    k = k + 1;
  END;

  FisherP=2*MIN(LEFTP, RIGHTP);

  Mid_P=2*MIN(midp_left, midp_Right);
  if mid_p>1 then MID_P=1;
  if mid_p<0 then MID_P=0;

END;

drop f lnfact z lnpi2 midp_left midp_right rightp p cell_a cell_b cell_c cell_d cell_r1 cell_r2 cell_c1 cell_c2 t pi pi2 e1 e2 e3 e4 e5
		loslop hislop  lnprob1 fisherp leftp leftp1 rightp1 sumcheck k n1 n2  ;*r1 r2 RR LL_RR UL_RR;
  ;
%MEND;

**********************************************EXAMPLE CODE********************************************************;
*MAKE SURE THAT THE ORDER OF THE CELLS (A,B,C,D) IS APPROPRIATE;

data x;
input AA BB CC DD;
cards;
-1 20 11 21
10 0 7 20
1 . 3 4
1 2 3 4
;
run;

DATA X;SET X;
%TWOBYTWO(A=AA,B=BB,C=CC,D=DD); 
RUN;
PROC PRINT;run;


