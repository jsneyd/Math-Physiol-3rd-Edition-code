
restart;
N:=1; # choose the number of barriers 
;
for j from 0 to N-1 do
eq||(j+1):= k||j*p||j-km||(j+1)*p||(j+1)-M;
od;
eq0:=koff*cibyKeq*x-koff*p0-M;
eq||(N+1):=k||N*p||N-koff*cebyKeq*x-M;
eq||(N+2):=x+sum('p||k',k=0..N)-1;
for j from 0 to N do
p||j:=solve(eq||j,p||j);
od;
M:=solve(eq||(N+1),M);
x:=solve(eq||(N+2),x);
M;
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

M;
factor((M));
# assume equal barrier heights
for j from 0 to N-1 do
 kpEmdG||j:=kpEmDG;
od;
for j from 1 to N do
kpEmdGm||j:=kpEmDG;
od;
M;
M:=factor(M);

dM:=denom(M);
factor(dM);
for j from 0 to 2*N do
dM||j:=factor(coeff(dM,EmVby||(2*N),j));
od;
NULL;
