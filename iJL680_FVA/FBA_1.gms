$INLINECOM /*  */
$set path output

OPTIONS
decimals = 8
solprint = on
reslim = 1000000
iterlim = 10000000
domlim = 0
limcol = 1000
limrow = 1000
optca = 0.0
optcr = 1E-9
work = 10000000
mip = cplex
;

SETS
i
/
$include "metabolites.txt"
/
j
/
$include "reactions.txt"
/

;

PARAMETERS
S(i,j)
/
$include "sij.txt"
/
LB(j)
/
$include "lower_bound.txt"
/
UB(j)
/
$include "upper_bound.txt"
/
;

VARIABLES
z
v(j)
R
;

POSITIVE VARIABLES
y(j)
m(j)
;

EQUATIONS
Obj
Stoich
con1
con2
con3
Biom1
Biom2
;

*Obj..            z =e= v('EX_cpd11416_norm2(e)')
Obj..            z =e= v('BIOMASS_Cl_DSM_WT_46p666M1')
;
Stoich(i)..      sum( j, S(i,j) * v(j) ) =e= 0
;
con1(j)..       v(j) =g= LB(j) + 0
;
con2(j)..       v(j) =l= UB(j) + 0
;

*v.fx('biomass_out') = 7.458;


MODEL MaxBiomass
/
Obj
Stoich
con1
con2
*Biom1
*Biom2
/
;

MaxBiomass.optfile = 1
;

SOLVE MaxBiomass USING LP MAXIMIZING z
;

file report /'%path%/FBARxnLevels.txt'/ ;
put report;

put 'min glucose level is:  ', z.l:0:10/;
put 'model stat is:  ',MaxBiomass.modelstat//;

alias(j,j1);
LOOP(j1,
    put j1.tl:0:60, "  ", v.l(j1):0:10/;
);
	



**** If Statement Example ****
file reportnonzero /'%path%/FBARxnLevels_NonZero1.txt'/ ;
file reportPositive /'%path%/FBARxnLevels_Positive1.txt'/ ;
file reportNegative /'%path%/FBARxnLevels_Negative1.txt'/ ;
file reportZero /'%path%/FBARxnLevels_Zero1.txt'/ ;
put reportnonzero;
alias(j,j2);
put reportnonzero;
LOOP(j2,
	   IF(v.l(j2)<>0,
			  put j2.tl:0:60,"  ",v.l(j2):0:10/;
			  );
	   );

LOOP(j2,
	   IF(v.l(j2)>0,
			  put reportPositive;
			  put j2.tl:0:60,"  ",v.l(j2):0:10/;
	   ELSEIF v.l(j2) lt 0,
			  put reportNegative;
			  put j2.tl:0:60,"  ",v.l(j2):0:10/;
	   ELSE
			  put reportzero;
			  put j2.tl:0:60,"  ",v.l(j2):0:10/;
			  );
	   );
