/**********************************************************************************************************************

									P-VALUE for SIR=1 and 95%CI OF SIR by mid-P exact method

**********************************************************************************************************************
Function: calculate 2 SIDED P-VALUE FOR 'SIR=1' or 'SIR=a nomimal value'

Developed by Minn M. Soe(Lead Epidemiologist, Surveillance Branch, NHSN, CDC), May 2025

validation: OpenEpi,PEPI

**********************************************************************************************************************

Input parameters:  OBS=observed value,  EXP=predicted value

Output: SIR=SIR, LLIMIT & ULIMIT=lower and upper limit of 95%CI of SIR

Note: No results will be generated if entry values are missing or impossible values.
Note: When SIR=0, the lower limit (LLIMIT) of 95%CI is set to missing by default. A user may change it to LLIMIT=0 if desired.

**********************************************************************************************************************
REFERENCE: 
1. JH ABRAMSON. WINPEPI PROGRAMS. DESCRIBE MANUAL (VERSION 2.42), PAGE-52. Available at 'http://www.brixtonhealth.com/pepi4windows.html'
2. Dean AG, Sullivan KM, Soe MM. OpenEpi: Open Source Epidemiologic Statistics for Public Health, Version. www.OpenEpi.com.
3. Rothman KJ, Boice JD Jr: Epidemiologic analysis with a programmable calculator. NIH Pub No. 79-1649. Bethesda, MD: National Institutes of Health, 1979;31-32.
4. Rothman KJ, Greenland S.  Modern Epidemiology, 2nd Edition.  Lippincott-Raven Publishers, Philadelphia, 1998.
5. GEOFFREY RC, SHU-YING Y. MID-P CONFIDENCE INTERVALS FOR THE POISSON EXPECTATION. STATISTICS IN MEDICINE, 13,2189-2203 (1994)
**********************************************************************************************************************/;


**********************************************************************************************************************

		**************************P-VALUE for SIR=1 ************************************

**********************************************************************************************************************;
%macro SIR(OBS=,EXP=);

e=2.718281828459045235; 

**********missing data handling and when &EXP<1;
IF  &OBS = . or &OBS<0  or mod(&OBS, 1) ^= 0 or &EXP = . or &EXP<1 THEN DO; 
midp=.;Byar_p=.;SIR=.;LLIMIT=.;ULIMIT=.;
end;
else do;

****when &OBServed numbers of events 'mu' are <=100;
if (&OBS<=100) /*AND (&EXP<=100)*/ THEN DO;

if (&OBS<&EXP) then do;
        k=&OBS-1; 
        do while (k>=0); 
    	numerator=(e**-&EXP)* (&EXP**k);
	
	    kk=k;
	    denom=1;
	    do while (kk>0);
			denom=denom*kk; kk=kk-1; *calculate the value of a factorial;
	    end;
	    subtotal=numerator/denom;	   
	    total=sum(total,subtotal);
    	k=k-1;
		end;

*first part of equation;
	num=0;  deno=0;  i=0;  aa=0;
	num=(e**-&EXP)* (&EXP**&OBS);
	deno=1;
	i=&OBS;
	do while (i>0); deno=deno*i; i=i-1;end; 
	
	aa=(num/deno)*0.5;
	
	*combine both parts of equation;
	MidP=2*(sum(aa,total));

end;


if (&OBS>=&EXP) then do;
        k=&OBS-1; 
        do while (k>=0); 
    	numerator=(e**-&EXP)* (&EXP**k);
	
	    kk=k;
	    denom=1;
	    do while (kk>0);
			denom=denom*kk; kk=kk-1; *calculate the value of a factorial;
	    end;
	    subtotal=numerator/denom;	   
	    total=sum(total,subtotal);
    	k=k-1;
		end;

*first part of equation;
	num=0;  deno=0;  i=0;  aa=0;
	num=(e**-&EXP)* (&EXP**&OBS);
	deno=1;
	i=&OBS;
	do while (i>0); deno=deno*i; i=i-1;end; 
	
	aa=(num/deno)*0.5;
	
	*combine both parts of equation;
	MidP=2*(1-(sum(aa+total)));

end;
if midp>1 then MidP=1;
if midp<0 then MidP=0;
END;


****when &OBServed numbers of events are >100, use Byar poisson approximation (Rothman KJ, Boice JD);
ELSE DO;

Byar_p=0; z=0; OBS_=0;
OBS_=&OBS;
if (&OBS<=&EXP) then OBS_=&OBS+1;
z=  sqrt(9*OBS_) * ( 1 - 1/(9*OBS_) - (&EXP/OBS_)**(1/3) );

if &OBS<=&EXP then Byar_p=2*probnorm(z);
else Byar_p=2*(1-probnorm(z));

if Byar_p>1 then Byar_p=1;
if Byar_p<0 then Byar_p=0;
midp=Byar_p;


END;


************************************************************************************************************

												95%CI

************************************************************************************************************;
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

*point estimate;
if &exp=>1 then do;
sir=&obs/&exp;

*Lower tail;
 v=0.5; dv=0.5; p=2.5/100;    	 	 
do while(dv>10**-5);
	dv=dv/2; 

	*component-A; 	
 	%fish(z=(1+&obs)*v/(1-v), x1=(&obs+1), x2=10**10); 

	*component-B;
	%fish2(z=(1+&obs)*v/(1-v), x1=&obs, x2=&obs); 

	if( (PoisP+PoisP2)>p) then v=v-dv; else v=v+dv;
end;

if &obs=0 then LLIMIT=.; *LEAVE LOWER LIMIT BLANK IF NO EVENT OCCURED;
else LLIMIT= ((1+&obs)*v/(1-v))/&exp;


*Upper tail;
 v=0.5; dv=0.5; p=2.5/100;    	 	 
do while(dv>10**-5);
	dv=dv/2; 

	*component-A; 	
 	%fish(z=(1+&obs)*v/(1-v), x1=0, x2=(&obs-1) ); 

	*component-B;
	%fish2(z=(1+&obs)*v/(1-v), x1=&obs, x2=&obs); 

	if( (PoisP+PoisP2)<p) then v=v-dv; else v=v+dv;
end;
ULIMIT= ((1+&obs)*v/(1-v))/&exp;

END;

END;

drop e Byar_p  k numerator kk denom subtotal total num deno i aa z OBS_ v dv p q tot s poisP poisP2 ;
%mend;


**********************************************EXAMPLE CODE********************************************************;
data x;
input obs exp;
cards;
1 2.3
0 3
;
run;

DATA X;SET X;
%sir(OBS=obs, exp=exp); 
RUN;
proc print;run;



