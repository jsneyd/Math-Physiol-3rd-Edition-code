from sympy import *

V0,V1,R = symbols("V0 V1 R")

V = V0 + V1*R
pp = 1 - (2*V/R)*(1 - exp(-R/(2*V)))*(1+V*R)
out = Poly(series(pp,R,0,3).removeO(),R)

cc = out.coeffs()
print(cc)
V0val = solve(cc[1],V0)[1]
print(V0val)

dum = cc[0].subs(V0,V0val)
print(solve(dum,V1))
