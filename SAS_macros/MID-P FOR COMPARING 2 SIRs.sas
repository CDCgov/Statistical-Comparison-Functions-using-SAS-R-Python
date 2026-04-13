/**********************************************************************************************************************

****************************Compare 2 standardized ratios (2 SIRs) based on mid-P exact method**********************
Function: compare 2 standardized ratios such as comparing 2 SIRs (or) similar poisson variates of standardized ratios such as SURs, SAARs, SMRs.
(based on the assumption that the distribution of two poisson variates conditional on their sum is binomial)

Developed by Minn M. Soe(Lead Epidemiologist, Surveillance Branch, NHSN, CDC), May 2025

Validation: Word Health Organization: Statistical Methods in Cancer Research

**********************************************************************************************************************

Data input layout:
O1= observed events in group1
E1= predicted events in group1
O2= observed events in group2
E2= predicted events in group2

Outputs:
midP= TWOTAILED MID-P value
RATIO= SIR in group2 / SIR in group1 = SIR2/SIR1 = (&O2/&E2)/(&O1/&E1)
LL= 95% lower bound of RATIO
UL= 95% Upper bound of RATIO

Note: When RATIO=0, the lower limit (LL) of 95%CI is set to missing by default. A user may change it to LL=0 if desired.

**********************************************************************************************************************

Reference 
1. Breslow NE, Day NE. Statistical methods in cancer research. Volume II--The design and analysis of cohort studies. IARC Sci Publ. 1987;(82):1-406. PMID: 3329634.
http://www.iarc.fr/en/publications/pdfs-online/stat/sp82/SP82_vol2-3.pdf
2. Erich L. Lehmann, Joseph P. Romano. Testing Statistical Hypotheses. Wiley, New York, 1959 (Third edition, 2005, Springer Texts in Statistics) 
3. John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
4. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com, updated 2013/03/21
**********************************************************************************************************************/;


%macro binom(O1=,E1=,O2=,E2=);

%macro BinP(x1=,x2=, PP=, N=);
q=&PP/(1-&PP);  k=0;  v = 1;  s=0;  tot=0;
    do while(k<=&N) ;
        tot=tot+v;
        if(k>=&x1 and k<=&x2) then s=s+v ;
        if(tot>10**30) then do;
			s=s/10**30; tot=tot/10**30; v=v/10**30;
		end;
        k=k+1; v=v*q*(&N+1-k)/k;
     end;
    binP=s/tot;
%mend;

%macro BinP2(x1=,x2=, PP=, N=);
q=&PP/(1-&PP);  k=0;  v = 1;  s=0;  tot=0;
    do while(k<=&N) ;
        tot=tot+v;
        if(k>=&x1 and k<=&x2) then s=s+v ;
        if(tot>10**30) then do;
			s=s/10**30; tot=tot/10**30; v=v/10**30;
		end;
        k=k+1; v=v*q*(&N+1-k)/k;
     end;
     binP2=(s/tot)*0.5;
%mend;


*input parameters;
IF (&O1=0 AND &O2=0) OR &O1=. OR &E1=. OR &O2=. OR &E2=. OR &E1<1 OR &E2<1  OR &O1<0 OR &O2<0 or mod(&O1, 1) ^= 0 or mod(&O2, 1) ^= 0 THEN DO;
midP=.; RATIO=.; LL=.; UL=.;
END;

ELSE DO;
vN=&O2+&O1; 
vP=&O2/vN; 


*******************MID-P VALUE FOR HYPOTHESIS TESTING;

RATIO1=&O1/&E1; RATIO2=&O2/&E2;
O=&O1+&O2;
T=&E1+&E2;

*take the larger of SIRs (TESTING FOR Ha: RR>1);
IF RATIO1>=RATIO2 THEN DO;
	p1=1- probbnml(&E1/T,O,&O1);
	p2=0.5* ( probbnml(&E1/T,O,&O1) - probbnml(&E1/T,O,&O1-1) );
	P3=P1+P2;
END;
ELSE DO;
	p1=1- probbnml(&E2/T,O,&O2);
	p2=0.5* ( probbnml(&E2/T,O,&O2) - probbnml(&E2/T,O,&O2-1) );
	P3=P1+P2;
END;

midP=2*p3;
IF midP>1 THEN midP=1;
IF midP<0 THEN midP=0;



*******************RATIO AND INTERVAL ESTIMATION (SIR2/SIR1, apriori assumption: SIR2>SIR1);
*Conditional maximum likelihood estimate of Ratio;
if &O1 ne 0 then do; *<-----------to prevent division by zero when ratio1=zero ;

RATIO=(&O2/&E2) / (&O1/&E1);

if (&O2=0) then DL = 0; 
else do;  p2=vP/2; vsL=0; vsH=vP;  p=2.5/100;
        do while((vsH-vsL)>10**-5);
			%BinP(x1=&O2+1, x2=vN, PP=P2, N=VN);
			%BinP2(x1=&O2,  x2=&O2, PP=P2, N=VN);
			 if (BinP+BinP2) >p  then do;vsH=p2; p2=(vsL+p2)/2 ;end;
			 else do; 			 vsL=p2; p2=(vsH+p2)/2 ;end;
		end;
        DL = P2;
end;     

if (&O2=vN) then DU = 1; 
else do;
        p3=(1+vP)/2; vsL=vP; vsH=1; p=2.5/100;
        do while((vsH-vsL)>10**-5);
			%BinP (x1=0, x2=&O2-1, PP=P3, N=VN);
			%BinP2(x1=&O2,  x2=&O2, PP=P3, N=VN);
			if (BinP+BinP2)  <p  then do; vsH=p3; p3=(vsL+p3)/2;end;
			else do; vsL=p3; p3=(vsH+p3)/2 ;end;
		end;
        DU = P3;
end;       

ll=(dl*&E1)/((1-dl)*&E2);
Ul=(DU*&E1)/((1-DU)*&E2);



***BECAUSE RATIO=sir2/ sir1,
WHEN sir2=0 then ratio OF 2 SIRS=0 AND ITS LOWER BOUND=.;
IF RATIO=0 THEN LL=.;

end;

END;

DROP vN  vP RATIO1 RATIO2 O T p1 p2 P3  DL vsL  vsH  p  q  k  v  s  tot  binP  binP2  DU;
run;

%mend;


**********************************************EXAMPLE CODE********************************************************;
data exe;
input obs2017 pred2017 obs2016 pred2016;
cards;
10 100 40 50
;
run;
data exe;set exe;
%binom(o1=obs2017, e1=pred2017, o2=obs2016, e2=pred2016);
run;
proc print;run;




