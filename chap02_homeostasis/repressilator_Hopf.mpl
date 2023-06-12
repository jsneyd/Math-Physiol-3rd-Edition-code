
restart;
eq1:=m/(1+y^4)-x;
eq2:=m/(1+z^4)-y;
eq3:=m/(1+x^4)-z;


with(linalg);
J:=array([[diff(eq1,x)-I*sqrt(w),diff(eq1,y),diff(eq1,z)],[diff(eq2,x),diff(eq2,y)-I*sqrt(w), diff(eq2,z)],[diff(eq3,x),diff(eq3,y), diff(eq3,z)-I*sqrt(w)]]);
dJ:=det(J);
dJ:=subs(y=x,z=x,dJ);
dJI:=simplify(coeff(dJ,I,1)/sqrt(w));
dJR:=factor(coeff(dJ,I,0));
dJI;
w:=solve(dJI,w);
simplify(dJR);
dJR:=numer(dJR);
z:=solve(eq3,z);
y:=solve(eq2,y);
factor(eq1);
factor(dJR);
m:=x^5+x;
factor(dJR);
x:=cx^(1/4);
factor(dJR);

evalf(solve(dJR,cx));
cx:=1;
m;
NULL;
