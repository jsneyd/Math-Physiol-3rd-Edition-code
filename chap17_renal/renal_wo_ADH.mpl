
restart; # the leading order problem
;
Cd:=(1+rd*Hd*(1-Qd)-DPd*Hd*y)/Qd;
Cs:=(P+DPd*Hd)*(1-y)/(Qd+Qa)-rd*Hd;
f:=Cd-Cs-DPd;
eq:= f *yp-rd;
y:=y0+rd*y1;
yp:=yp0+rd*yp1; 
  
Eq:=convert(series(eq,rd),polynom):
Eq0:=coeff(Eq,rd,0);
Y0:=solve(Eq0,y0);
simplify(subs(Qd=-Qa,Y0));
simplify(subs(Qd=1,Y0));
Qa:=solve(%,Qa); simplify(1+Qa);
simplify(Y0-1);
P:=8/10; DPd := 3/100; Hd:=1/10;  
Y0;
Qa;  # be sure Qa<0;
plot(Y0,Qd = -Qa..1);
subs(Qd=-Qa,Y0)
;
(DPd*Hd + P)/(-DPd + 1);
restart;  # higher order corrections
;
Cd:=(1+rd*Hd*(1-Qd)-DPd*Hd*y)/Qd;
Cs:=(P+DPd*Hd)*(1-y)/(Qd+Qa)-rd*Hd;
f:=Cd-Cs-DPd;
eq:= f/(Qa+1) *yp-rd;
Qd:=(x-1)*Qa+x;
subs(x=1,Qd); subs(x=0,Qd);  diff(Qd,x)
;
# powers series;
y:=y0+rd*y1; #+rd^2*y2;
yp:=yp0+rd*yp1;#+rd^2*yp2;
Qa:=Qa0+Qa1*rd;#+Qa2*rd^2;
Eq:=convert(series(eq,rd),polynom):
Eq0:=coeff(Eq,rd,0);
Y0:=solve(Eq0,y0);
simplify(subs(x=0,Y0));
simplify(subs(x=1,Y0));
Qa0:=solve(%,Qa0);
simplify(Y0);
simplify(subs(x=1,Y0));
plot(subs(P=8/10,DPd = 5/100,Hd=1/10,Y0),x=0..1);
# next order
Eq1:=coeff(Eq,rd,1):
y0:=Y0;
yp0:=diff(y0,x):
Y1:=solve(Eq1,y1): y1:=Y1;
subs(x=0,Y1);
subs(x=1,Y1);
Qa1:=solve(subs(x=1,Y1),Qa1);
plot([subs(P=8/10,DPd = 5/100,Hd=1/10,rd=1/100,y),subs(P=8/10,DPd = 5/100,Hd=1/10,rd=1/100,Qd),x=0..1]);
 plot([subs(P=8/10,DPd = 5/100,Hd=1/10,rd=1/100,y),subs(P=8/10,DPd = 5/100,Hd=1/10,rd=1/100,Cs),x=0..1]);

NULL;
NULL;
NULL;
