
NULL;
restart;
eq1:=-a1*p00+km1*p10+k2*p01;
eq2:=k0*ci*p00-a2*p10+km2*p01+k2*p11;
eq3:=km3*ce*p00+k1*p10-a3*p01+km1*p11;
eq4:=km3*ce*p10+k0*ci*p01-a4*p11;
a1:=k0*ci+km3*ce;
a2:=km1+k1+km3*ce;
a3:=k2+km2+k0*ci;
a4:=k2+km1;
simplify(eq1+eq2+eq3+eq4);
eq5:=p00+p10+p01+p11-1;
solve({eq1,eq2,eq3,eq5},{p00,p01,p10,p11})
;
assign(%);
J:=k1*p10-km2*p01;
factor(J);
limit(J,ci=infinity);

