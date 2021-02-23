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
;

EQUATIONS
Obj
Stoich
con1
con2
con3
Biom1

;

Obj..            z =e= sum( j , c(j)*v(j) )
;
Stoich(i)..      sum( j, S(i,j) * v(j) ) =e= 0
;
con1(j)..       v(j) =g= LB(j)
*- y(j) 
;
con2(j)..       v(j) =l= UB(j)
*+ m(j)
;
*con3..       sum(j, y(j)) + sum(j,m(j)) =e= R
*;
Biom1..      v('Biomass_c') =g= 0.241561 -  0.241561 * 0.001;
*Biom2..      v('biomass_out') =g= 8.97335203471916 -  8.97335203471916 * 0.525;



MODEL MaxBiomass
/
Obj
Stoich
con1
con2
Biom1
*con3

/
;

MaxBiomass.optfile = 1
;

file report /"%path%/fva_aero_solcheck1.txt"/;
file block /"%path%/fva_blocked1.txt"/;
file invar /"%path%/fva_invar1.txt"/;
file var /"%path%/fva_var1.txt"/;
put report

alias(j,j2);
loop(j2,
  c(j)=0;
  c(j2)=1;
 SOLVE MaxBiomass USING LP MAXIMIZING z;
  model_stat = MaxBiomass.modelstat;
  PUT j2.tl:40, "   ", z.l:10:5,"    ",model_stat/;
  high(j2) = z.l;

  SOLVE MaxBiomass USING LP MINIMIZING z;
  model_stat = MaxBiomass.modelstat;
 PUT j2.tl:40, "   ",z.l:10:5,"   ", model_stat/;
  low(j2) = z.l; 

);
putclose report
put block
PUT "Blocked reactions"/;
loop(j$(high(j) lt epsilon and low(j) gt -epsilon),
        put j.tl:30:0, high(j):0:8,"  ",low(j):0:8/;
);
putclose block

put invar
put /"Invariant reactions"/;
loop(j$((high(j) gt epsilon or low(j) lt -epsilon) and ((high(j) - low(j)) lt epsilon)),
        put j.tl:30:0, high(j):0:8,"  ",low(j):0:8/;
);
putclose invar

put var
put /"Variable-range reactions"/;
loop(j$((high(j) - low(j)) gt epsilon),
        put j.tl:30:0, high(j):0:8,"  ",low(j):0:8/;
);
putclose var

putclose report
