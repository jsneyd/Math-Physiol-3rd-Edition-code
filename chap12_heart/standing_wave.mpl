
restart;
NULL;
NULL;
ph:=a1*exp(l*x) + a2*exp(-l*x);
phii:=ratc*ph+b;  # ratc=rc/(rc+re);
phe:=-(1-ratc)*ph+b;
NULL;
eq1:=subs(x=L,diff(phe,x))-m*subs(x=0,diff(phe,x));
eq2:=subs(x=L,diff(phii,x))-m*subs(x=0,diff(phii,x));
eq3:=subs(x=L, phe)-m*subs(x=0, phe);
eq4:=subs(x=L, phii)-m*subs(x=0, phii)+rgbyrc*subs(x=L,diff(phii,x));
solve({eq1,eq2,eq4},{a1,a2,b});
assign(%);
a2:=m-exp(l*L);
a1;
b;
L:=log(E)/l;
b;
simplify(b);
eq3;
lam:=solve(eq3,l);
simplify(subs(l=lam,b));
simplify(a1);
a2;
#now find the standing wave solution
;
Vi0:=A*phii;
Ve0:=A*phe;
Vim1:=B*m*subs(m=1/n,phii)+1+C;
Vem1:=B*m*subs(m=1/n,phe)+C;
# boundary conditions;
eq5:=subs(x=L,Vem1)-subs(x=0,Ve0); #Ve continuous
;
eq6:=subs(x=L,diff(Vem1,x))-subs(x=0,diff(Ve0,x)); # continuous currentdir
;
eq7:=subs(x=L,diff(Vim1,x))-subs(x=0,diff(Vi0,x));  # continuous current 
;
eq8:=subs(x=L,Vim1)-subs(x=0,Vi0)+rgbyrc*subs(x=L,diff(Vim1,x)); # jump condition 
;
solve({eq5,eq6,eq8},{A,B,C});
assign(%);
n:=m;
A; A/B;
C;
simplify(eq7);
simplify(subs(x=0,1/(Vi0-Ve0))-2);
simplify(subs(x=0,(Vi0-Ve0))+subs(x=L,Vim1-Vem1));
subs(l=LAM/(ratc*rgbyrc),A);
subs(l=LAM/(ratc*rgbyrc),C);
NULL;
