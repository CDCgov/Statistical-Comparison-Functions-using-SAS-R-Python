/************************************************************************************************************************

								95%CI of Incidence Density Rate by mid-P exact method

*************************************************************************************************************************
Function: calculate 95%CI of AN INCIDENCE DENSITY RATE

Developed by Minn M. Soe(Lead Epidemiologist, Surveillance Branch, NHSN, CDC), May 2025

validation: OpenEpi,PEPI

***********************************************************************************************************************

Input parameters:  numer=numerator,  denom=denominator.

output: incidence density rate per 1000 person-time & lower and upper limit of 95%CI: rate_l & rate_u.

Note: When rate=0, the lower limit (rate_l) of 95%CI is set to missing by default. A user may change it to rate_l=0 if desired.


***********************************************************************************************************************

Reference 
1. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version 2.3.1. www.OpenEpi.com
2. Rothman KJ,  Boice JD  Jr:  Epidemiologic analysis with a programmable calculator.   NIH Pub No. 79-1649.  Bethesda, MD:  National Institutes of Health, 1979;31-32.
3. Rothman KJ, Greenland S.  Modern Epidemiology, 2nd Edition.  Lippincott-Raven Publishers, Philadelphia, 1998.

***********************************************************************************************************************/;


%macro rateCIComp(numer=,denom=);
*input parameters;
****HANDLING OF MISSING/IMPOSSIBLE VALUES AND RELATION;
IF (&numer=0 and &denom=0)  or &numer=. or &numer<0 or &denom=. or &denom<=0 or &numer>&denom  or mod(&numer, 1) ^= 0 or mod(&denom, 1) ^= 0 THEN DO;
rate=.;rate_l=.; rate_u=.;
END;

ELSE DO;
rate=(&numer/&denom)*1000;

%macro fish(z=,x1=,x2=);
q=1;  tot=0;  s=0;  k=0;
    do while(k<&z or q>(tot*10**-10)); 
        tot=tot+q;
        if(k>=&x1 and k<=&x2) then s=s+q ;
       	if(tot>10**30)then do;
			s=s/10**30; tot=tot/10**30; q=q/10**30;
		end;
        k=k+1; q=q*&Z/k;
     end;
  	 poisP=s/tot;
%mend;
%macro fish2(z=,x1=,x2=);
q=1;  tot=0;  s=0;  k=0;
    do while(k<&z or q>(tot*10**-10)); 
        tot=tot+q;
        if(k>=&x1 and k<=&x2) then s=s+q ;
       	if(tot>10**30)then do;
			s=s/10**30; tot=tot/10**30; q=q/10**30;
		end;
        k=k+1; q=q*&Z/k;
     end;
  	 poisP2=0.5*(s/tot);
%mend;

*Lower tail;
	if &numer=0 then rate_l=.;*new update;
	else do;
	v=0.5; dv=0.5; p=2.5/100;    	 	 
	do while(dv>10**-5);
		dv=dv/2; 
	 	%fish(z=(1+&numer)*v/(1-v), x1=(&numer+1), x2=10**10); 
		%fish2(z=(1+&numer)*v/(1-v), x1=&numer, x2=&numer); 
		if( (PoisP+PoisP2)>p) then v=v-dv; else v=v+dv;
	end;
	if &denom > 0 then rate_l=(((1+&numer)*v/(1-v))/&denom)*1000;
	end;


*Upper tail;
	v=0.5; dv=0.5; p=2.5/100;    	 	 
	do while(dv>10**-5);
		dv=dv/2; 
	 	%fish(z=(1+&numer)*v/(1-v), x1=0, x2=(&numer-1) ); 
		%fish2(z=(1+&numer)*v/(1-v), x1=&numer, x2=&numer); 
		if( (PoisP+PoisP2)<p) then v=v-dv; else v=v+dv;
	end;
	if &denom > 0 then rate_u=(((1+&numer)*v/(1-v))/&denom)*1000;
	drop v dv p q tot k s poisP poisP2;
END;
%mend;


***************************************************example code*****************************************************;
data z;
input obs persontime;
cards;
0 10
10 10
;
run;
DATA z;SET z;
%rateCIcomp(numer=obs, denom=persontime); 
RUN;
proc print;run;
