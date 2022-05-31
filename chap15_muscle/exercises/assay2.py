from sympy import *

A,g2 = symbols("A g2")

vl = g2/2 + A/24;
vs = g2/sqrt(2) - g2*g2/(2*A);

dd = vl - vs
g2m = solve(diff(dd,g2),g2)
print(simplify(dd.subs(g2,g2m[0])))
