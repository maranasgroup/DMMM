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
epsilon /1e-6/
model_stat
c(j)
high(j)
low(j)
;

VARIABLES
z
v(j)
R
;

POSITIVE VARIABLES
y(j)
m(j)
n(j)
;

EQUATIONS
Obj
Stoich
con1
con2
con3
con4

;

Obj..            z =e= sum( j , m(j) ) + sum( j , n(j) )
;
Stoich(i)..      sum( j, S(i,j) * v(j) ) =e= 0
;
con1(j)..       v(j) =g= LB(j) - m(j)

;
con2(j)..       v(j) =l= UB(j) + n(j)

;
m.fx('EXCH_accoa(e)') = 0;
n.fx('EXCH_accoa(e)') = 0;
m.fx('EXCH_ser(e)') = 0;
n.fx('EXCH_ser(e)') = 0;
m.fx('EXCH_gly-l(e)') = 0;
n.fx('EXCH_gly-l(e)') = 0;
m.fx('EXCH_asp(e)') = 0;
n.fx('EXCH_asp(e)') = 0;
m.fx('EXCH_met-l(e)') = 0;
n.fx('EXCH_met-l(e)') = 0;
m.fx('EXCH_phe-l(e)') = 0;
n.fx('EXCH_phe-l(e)') = 0;
m.fx('EXCH_cit(e)') = 0;
n.fx('EXCH_cit(e)') = 0;
m.fx('EXCH_fum(e)') = 0;
n.fx('EXCH_fum(e)') = 0;
m.fx('EXCH_pep(e)') = 0;
n.fx('EXCH_pep(e)') = 0;


MODEL MaxBiomass
/
Obj
Stoich
con1
con2


/
;

MaxBiomass.optfile = 1
;

SOLVE MaxBiomass USING LP MINIMIZING z
;

file report /'%path%/FBARxnLevels.txt'/ ;
put report;

put 'min glucose level is:  ', z.l:0:10/;
put 'model stat is:  ',MaxBiomass.modelstat//;

alias(j,j1);
LOOP(j1,
    put j1.tl:0:60, "  ", v.l(j1):0:10,"    +",n.l(j1):0:10,"     -",m.l(j1):0:10/;
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
