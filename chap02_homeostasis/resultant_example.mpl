
restart;
eq1:=a*y-1/y+z/x;  #using x instead of mu
;
eq2:=a*y+1/y+1/x-2;
R1:=factor(resultant(numer(eq1),numer(eq2),y));
NULL;
R2:=factor(resultant(numer(eq1),numer(eq2),x));
NULL;
