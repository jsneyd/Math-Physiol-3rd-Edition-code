
restart;
# step one: Find J as a function of x
NULL;
#case n=3 
 
eq1 := k0*c0*x - km1*c1 - k1*c1 + km2*c2;
eq2 := k1*c1 - km2*c2 - k2*c2 + km3*c3*x;
solve({eq1,eq2},{c1,c2}):
assign(%);
J := simplify(k0*c0*x - km1*c1);
NULL;
restart;
# Find J as a function of x
NULL;
#case n=4; 
eq1 := k0*c0*x - km1*c1 -  k1*c1 + km2*c2;
eq2 := k1*c1 - km2*c2 - k2*c2 + km3*c3;
eq3 := k2*c2 - km3*c3 - k3*c3 + km4*c4*x;
solve({eq1,eq2,eq3},{c1,c2,c3}):
assign(%);
J := simplify(k0*c0*x - km1*c1);
#step 2: eliminate x
;
restart;
#case n=3;
eq1 := k0*c0*x - km1*c1 - k1*c1 + km2*c2;
eq2 := k1*c1 - km2*c2 - k2*c2 + km3*c3*x;
# add conservation
cons := x+c1+c2-1;
solve({eq1,eq2,cons},{c1,c2,x}):
assign(%);
J := simplify(k0*c0*x - km1*c1);
NULL;
restart;
#case n=4;
eq1 := k0*c0*x - km1*c1 -  k1*c1 + km2*c2;
eq2 := k1*c1 - km2*c2 - k2*c2 + km3*c3;
eq3 := k2*c2 - km3*c3 - k3*c3 + km4*c4*x;
cons :=x+c1+c2+c3-1;
solve({eq1,eq2,eq3,cons},{c1,c2,c3,x}):
assign(%);
J := simplify(k0*c0*x - km1*c1);
NULL;
NULL;
