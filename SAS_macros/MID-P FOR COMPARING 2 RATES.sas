
/**********************************************************************************************************************;

**********************************COMPARE 2 INCIDENCE DENSITY RATES by mid-P exact method***********************************

Function: compare 2 incidence density rates (based on the fact that the distribution of two poisson variates conditional on their sum is binomial)

Developed by Minn M. Soe(Lead Epidemiologist, Surveillance Branch, NHSN, CDC), May 2025

Validation: WinPEPI, OpenEpi, Propective study of dietary fat and risk of prostate cancer. J national cancer institute. 1993, 85:1571-79.

**********************************************************************************************************************

Data input layout:
O1= observed events in group-1
PT1=person-time in group-1
O2= observed events in group-2
PT2= person-time in group-2

OUTPUT:
mid-P=significance test for comparing between 2 rates
RATE_RATIO = Rate in group-2 / Rate in group-1
LL=lower limit of 95%CI of rate ratio
UL=Upper limit of 95%CI of rate ratio

Note: When RATE_RATIO=0, the lower limit (LL) of 95%CI is set to missing by default. A user may change it to LL=0 if desired.

**********************************************************************************************************************

Reference 
1. Erich L. Lehmann, Joseph P. Romano. Testing Statistical Hypotheses. Wiley, New York, 1959(Third edition, 2005, Springer Texts in Statistics) 
2. John Pezzullo. Interactive Statistical Calculation Pages. www.Statpages.com.
3. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com, updated 2013/03/21
4. Austin, Harland (Emory). Epidemiology method-II: Statistical issues for Density type follow-up studies.


**********************************************************************************************************************/;

%macro TWORATES(O1=,PT1=,O2=,PT2=);

*STAT FUNCTION FOR INTERVAL CALCULATION;
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
****HANDLING OF MISSING/ EXTREME VALUES;
IF (&O1=0 and &O2=0) or &O1=. or &O1<0 or &PT1=. or &PT1<=0 or &O2=. or &O2<0 or &PT2=. or &PT2<=0 or &O1>&PT1 or &O2>&PT2 or mod(&O1, 1) ^= 0 or mod(&O2, 1) ^= 0 or mod(&PT1, 1) ^= 0 or mod(&PT2, 1) ^= 0 THEN DO;
MID_P=.;RATE_RATIO=.; LL=.;UL=.;
END;

ELSE DO;
vN=&O2+&O1; 
vP=&O2/vN; 


****************MID-P FOR HYPOTHESIS TESTING;
RATIO1=&O1/&PT1; RATIO2=&O2/&PT2;
O=&O1+&O2;
T=&PT1+&PT2;

*take the larger of RISKs (TESTING FOR Ha: RR>1);
IF RATIO1>=RATIO2 THEN DO;
	p1=1- probbnml(&PT1/T,O,&O1);
	p2=0.5* ( probbnml(&PT1/T,O,&O1) - probbnml(&PT1/T,O,&O1-1) );
	P3=P1+P2;
END;
ELSE DO;
	p1=1- probbnml(&PT2/T,O,&O2);
	p2=0.5* ( probbnml(&PT2/T,O,&O2) - probbnml(&PT2/T,O,&O2-1) );
	P3=P1+P2;
END;

MID_P=2*p3;
IF MID_P>1 THEN MID_P=1;
IF MID_P<0 THEN MID_P=0;


*****************RATE_RATIO AND INTERVAL ESTIMATION;
*Conditional maximum likelihood estimate of Ratio= rate2/rate1;
IF &O1 NE 0 THEN DO; *PREVENT DIVISION BY ZERO WHEN RATE RATIO1=0;

RATE_RATIO=(&O2/&PT2) / (&O1/&PT1);


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

ll=(dl*&PT1)/((1-dl)*&PT2);
Ul=(DU*&PT1)/((1-DU)*&PT2);


***BECAUSE RATE_RATIO=RISK2/ RISK1....
WHEN RISK GROUP2=0 then rate-ratio=0;
IF RATE_RATIO=0 THEN LL=.;

end;
END;

drop O  P3  RATIO1 RATIO2 T binP binP2 k p  p1 p2 q s tot v vN vP vsH vsL DL DU;
%mend;


**********************************************EXAMPLE CODE********************************************************;
data A;
input Obs1 PTime1 Obs2 PTime2;
cards;
10 200 11 210
;
run;
DATA A;SET A;
%TWORATES(O1=OBS1, PT1=PTIME1, O2=OBS2, PT2=PTIME2); 
RUN;
proc print;run;


                               

