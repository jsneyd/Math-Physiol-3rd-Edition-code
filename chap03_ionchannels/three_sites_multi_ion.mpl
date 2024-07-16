
NULL;
restart; # there are 3 binding sites and 4 barriers, multi-ion
;
# write out the equations
eq1:=-a1*p1+km4ce*p2+k0ci*p4+k2*p6+km2*p8;
eq2:=k3*p1-a2*p2+k0ci*p3+km2*p7;
eq3:=km1*p2-a3*p3+k3*p4;
eq4:=km1*p1+km4ce*p3-a4*p4+k2*p7;
eq5:=-a5*p5+km4ce*p6+k0ci*p8;
eq6:=km3*p1+k3*p5-a6*p6+k0ci*p7;
eq7:=k1*p2+km1*p4+km1*p6-a7*p7+k3*p8;
eq8:=k1*p1+km1*p5+km4ce*p7-a8*p8;
a1:= k1+ km1+k3+km3;
a2:=k1 + km1 + km4ce;
a3:=km4ce+k0ci;
a4:=k0ci+k3+km1;
a5:=k3+km1;
a6:=km4ce+km1+k2;
a7:=km4ce+k0ci+k2+km2;
a8:=km2+k0ci+k3;   
 eq9:=p1+p2+p3+p4+p5+p6+p7+p8-1;
 simplify(eq1+eq2+eq3+eq4+eq5+eq6+eq7+eq8);
N:=4; # number of barriers
;
solve({eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq9},{p1,p2,p3,p4,p5,p6,p7,p8}):
assign(%);
J:=k0ci*(p3+p4+p7+p8)-k1*(p1+p2+p5+p6):

J2:=km4c3*(p2+p3+p6+p7) - k3*(p1+p5+p4+p8):
J3:=k2*(p6+p7)-km3*(p1+p4):
J4:=k1*(p2+p1)-km2*(p8+p7):
km3:=km1;;
km2:=km1;
k3:=k1;
k2:=k1;
factor(numer(J3));
simplify(J3-J4);
