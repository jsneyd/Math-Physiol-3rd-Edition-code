
NULL;
NULL;
NULL;
;
NULL;
NULL;
restart;
N:=5;
NULL;

NULL;
eq0 := k0*ci - km1*c1 - M;
NULL;
# 
# 
eqNN := k||(N)*c||N - km||(N+1)*ce - M;
NULL;
;
for i from 1 to N-1 do 
  eq(i):= k||(i)*c||i - km||(i+1)*c||(i+1) - M;
od;

c1:=solve(eq0,c1);

for j from 2 to N do
c||j:=solve(eq(j-1),c||j);
od;

NULL;
M:=solve(eqNN,M);

NULL;
