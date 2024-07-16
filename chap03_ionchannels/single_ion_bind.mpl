
restart;# this is code to calculate the flux through a single-ion binding ion channel model
;
N:=2; # choose the number of binding sites =N+1 barriers 
;
for j from 1 to N-1 do
eq||(j):= k||j*p||j-km||(j+1)*p||(j+1)-J;
od;
eq0:=k0*ci*x-km1*p1-J;
eq||(N):=k||(N)*p||(N)-km||(N+1)*ce*x-J;
eq||(N+1):=x+sum('p||k',k=1..N)-1;
for j from 1 to N do
p||j:=solve(eq||(j-1),p||j);
od;
J:=solve(eq||(N),J);
x:=solve(eq||(N+1),x);
J; simplify(limit(J,ci=infinity));simplify(limit(J,ci=0));
#now add dG and V dependence
;
for j from 0 to N do
k||j:=kbar||j/EmVby||(2*(N+1));
od;
for j from 1 to N+1  do
km||j:=kbarm||(j)*EmVby||(2*(N+1));
od;
J:=factor(J);;
limit(J,ci=infinity);

