
restart;
N:=1; # choose the number of barriers 
;
for j from 0 to N-1 do
eq||(j+1):= k||j*p||j-km||(j+1)*p||(j+1)-J;
od;
eq0:=koff*cibyKeq*x-koff*p0-J;
eq||(N+1):=k||N*p||N-koff*cebyKeq*x-J;
eq||(N+2):=x+sum('p||k',k=0..N)-1;
for j from 0 to N do
p||j:=solve(eq||j,p||j);
od;
J:=solve(eq||(N+1),J);
x:=solve(eq||(N+2),x);
J;
#now add dG and V dependence
;
for j from 0 to N-1 do
k||j:=kpEmdG||j/EmVby||(2*N);
od;
for j from 1 to N do
km||j:=kpEmdGm||j*EmVby||(2*N);
od;
km0:=koff;
k||(N):=koff;

J;
factor(J);
# assume equal barrier heights
for j from 0 to N-1 do
 kpEmdG||j:=kpEmDG;
od;
for j from 1 to N do
kpEmdGm||j:=kpEmDG;
od;
J;
J:=factor(J);

dJ:=denom(J);
factor(dJ);
for j from 0 to 2*N do
dJ||j:=factor(coeff(dJ,EmVby||(2*N),j));
od;
limit(J,cibyKeq=infinity)
;

limit(J,cebyKeq=infinity);
NULL;
NULL;
NULL;
NULL;
NULL;
NULL;
NULL;
NULL;
NULL;
