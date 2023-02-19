
restart;  # multi-ion single barrier
;
# could use matrices here as well
eq1:=-2*koff*p1+koff*cebyKeq*p2+koff*cibyKeq*p4;
eq2:=koff*p1-(koff*cebyKeq+koff+k42)*p2+koff*cibyKeq*p3+k24*p4;
eq3:=koff*p2-(koff*cibyKeq+koff*cebyKeq)*p3+koff*p4;
eq4:=p1+p2+p3+p4-1;
solve({eq1,eq2,eq3,eq4},{p1,p2,p3,p4});
assign(%);
J:=p2*k42-p4*k24;
factor(J);
k42:=kpEmdG/EmVby2;
k24:=kpEmdG*EmVby2;
numer(J);
dJ:=denom(J);
for j from 1 to 3 do
factor(coeff(dJ,EmVby2,j));
od
;
Js:=subs(EmVby2=1,J);
factor(Js);
simplify(limit(J,cibyKeq=infinity));
NULL;
