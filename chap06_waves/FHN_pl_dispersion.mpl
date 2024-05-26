
restart;
w2:=A1*exp(l1*x)+A2*exp(l2*x)+A3*exp(l3*x);
w1:=1+B1*exp(l1*x)+B2*exp(l2*x)+B3*exp(l3*x);
v2:=-c*diff(w2,x);
v1:=-c*diff(w1,x);
eq1:=subs(x=x1,w2-w1);
eq2:=subs(x=x1,v2-v1);
eq3:=subs(x=x1,(diff(v2-v1,x)));
eq4:=subs(x=0,w1)-subs(x=x2,w2);
eq5:=subs(x=0,v1)-subs(x=x2,v2);
eq6:=subs(x=0,diff(v1,x))-subs(x=x2,diff(v2,x));
eq7:=subs(x=0,v1)-alp;
eq8:=subs(x=x1,v1)-alp;
solve({eq1,eq2,eq3,eq4,eq5,eq6},{A1,A2,A3,B1,B2,B3});
assign(%);
eq7;
eq8;
