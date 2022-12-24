
restart;
dex1:=km1*ni3*x2+k4*y1-(k1*ci+km4)*x1;
dex2 := km2*y2 + ci*k1*x1 - (km1*ni3 + k2)*x2;
dey1:=km4*x1+k3*ne3*y2-(k4+km3*ce)*y1;
eq4:=x1+x2+y1+y2-1;
solve({dex1,dex2,dey1,eq4},{x1,x2,y1,y2});
assign(%);
J:=k4*y1-km4*x1;
factor(numer(J));
NULL;
NULL;
NULL;
