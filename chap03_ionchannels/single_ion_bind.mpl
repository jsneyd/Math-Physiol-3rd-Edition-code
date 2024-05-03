
restart;# this is code to calculate the flux through a single-ion binding ion channel model
;
N:=1; # choose the number of barriers =N+1 binding sites
;
for j from 1 to N do
eq||(j):= k||j*p||j-km||(j+1)*p||(j+1)-J;
od;
eq0:=k0*ci*x-km1*p1-J;
eq||(N+1):=k||(N+1)*p||(N+1)-km||(N+2)*ce*x-J;
eq||(N+2):=x+sum('p||k',k=1..N+1)-1;
for j from 1 to N+1 do
p||j:=solve(eq||(j-1),p||j);
od;
J:=solve(eq||(N+1),J);
x:=solve(eq||(N+2),x);
J; limit(J,ci=infinity);
#now add dG and V dependence
;
for j from 0 to N do
k||j:=kpEmdG||j/EmVby||(2*(N+1));
od;
for j from 1 to N+1 do
km||j:=kpEmdGm||j*EmVby||(2*(N+1));
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
limit(J,ci=infinity)
;

limit(J,ceb=infinity);

