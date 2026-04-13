/**********************************************************************************************************************

								95%CI of single proportion by mid-P exact method 

*************************************************************************************************************************

Function: calculate 95%CI FOR SINGLE PROPORTION

Developed by Minn M. Soe(Lead Epidemiologist, Surveillance Branch, NHSN, CDC), May 2025

validation: OpenEpi,PEPI

**********************************************************************************************************************
Input parameters:  vx=numerator,  vn=denominator.

Output: vp=proportion, dl=lower limit of 95%CI, du=upper limit.

Note: When vp=0, the lower limit (dl) of 95%CI is set to missing by default. A user may change it to dl=0 if desired.

**********************************************************************************************************************

Reference 
1. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version 2.3.1. www.OpenEpi.com
2. Rothman KJ,  Boice JD  Jr:  Epidemiologic analysis with a programmable calculator. NIH Pub No. 79-1649.  Bethesda, MD:  National Institutes of Health, 1979;31-32.
3. Fleiss JL.  Statistical Methods for Rates and Proportions, 2nd Ed. John Wiley & Sons, New York, 1981.
4. Newcombe RG.  Two-sided confidence intervals for the single proportion: comparison of seven methods.  Statistics in Medicine 1988;17:857-872.
5. Stein Vollset. CONFIDENCE INTERVALS FOR A BINOMIAL PROPORTION STATISTICS IN MEDICINE,12,809-824 (1993)*/

**********************************************************************************************************************;

%macro binom(vX=,vN=);
*handing missing values/ illogical values;*<--------new update;
if &vX=. or &vN in (.  0) or &vX<0 or &vN<0 or  &vN<&vx  or mod(&vX, 1) ^= 0 or mod(&vN, 1) ^= 0 then do;
vp=.;dl=.;du=.;
end;
else do;
%macro BinP(x1=,x2=);
q=p2/(1-p2);  k=0;  v = 1;  s=0;  tot=0;
    do while(k<=&vN) ;
        tot=tot+v;
        if(k>=&x1 and k<=&x2) then s=s+v ;
        if(tot>10**30) then do;
			s=s/10**30; tot=tot/10**30; v=v/10**30;
		end;
        k=k+1; v=v*q*(&vN+1-k)/k;
     end;
    binP=s/tot;
%mend;
%macro BinP2(x1=,x2=);
q=p2/(1-p2);  k=0;  v = 1;  s=0;  tot=0;
    do while(k<=&vN) ;
        tot=tot+v;
        if(k>=&x1 and k<=&x2) then s=s+v ;
        if(tot>10**30) then do;
			s=s/10**30; tot=tot/10**30; v=v/10**30;
		end;
        k=k+1; v=v*q*(&vN+1-k)/k;
     end;
    binP2=(s/tot)*0.5;
%mend;


vP=&vX/&vN;

if (&vX=0) then DL = .; *<--------new update;
else do;
        p2=vP/2; vsL=0; vsH=vP;  p=2.5/100;
        do while((vsH-vsL)>10**-5);
			%BinP(x1=&vX+1, x2=&vN);
			%BinP2(x1=&vX, x2=&vX);
			 if (BinP+BinP2) >p  then do; vsH=p2; p2=(vsL+p2)/2 ;end;
			else do; vsL=p2; p2=(vsH+p2)/2 ;end;
		end;
        DL = p2;
end;       
   

if (&vX=&vN) then DU = 1; 
else do;
        p2=(1+vP)/2; vsL=vP; vsH=1; p=2.5/100;
        do while((vsH-vsL)>10**-5);
			%BinP(x1=0, x2=&vX-1);
			%BinP2(x1=&vX, x2=&vX);
			if (BinP+BinP2) <p then do; vsH=p2; p2=(vsL+p2)/2;end;
			else do; vsL=p2; p2=(vsH+p2)/2 ;end;
		end;
        DU = p2;
end; 

drop p2 vsl vsh p q k v s tot binp binp2; 
END;
%mend;


*--------------------------------example data set---------------------------------------;
data test;
input numerator denominator;
cards;
5 10
;
run;
data test;set test;
%binom(vx=numerator, vn=denominator);
run;
proc print; run;


