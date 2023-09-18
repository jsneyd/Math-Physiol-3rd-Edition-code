
#Baroreceptor loop
;
restart;
eq1:=F*Cd*Pv-Q;
eq2:=Pa-Pv-Q*Rs;
eq3:=Ca*Pa-Va;
eq4:=Cv*Pv-Vv;
eq5:=Va+Vv-Vt;
eq6:=(O2a-O2v)*Q-M;
Rs:=Kp*Smx/(Kp+Pa)+A*O2v;;
solve({eq1,eq2,eq3,eq4,eq5,eq6},{Q, F, Va, Vv, Pv, O2v});
assign(%);
factor(numer(F));
factor(numer(Q));
F/Q;

Dpa:=diff(Q,Pa);
NDpa:=numer(Dpa);

factor(coeff(NDpa,Smx,1));
factor(numer(Q));
NULL;
