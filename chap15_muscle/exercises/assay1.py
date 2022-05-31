from sympy import *

V0,V1,G = symbols("V0 V1 G")

V = V0 + V1*G
pp = 1 - (2*V*G)*(1)*(1+V/G)
out = Poly(series(pp,G,0,2).removeO(),G)

cc = out.coeffs()
print(cc)
V0val = solve(cc[1],V0)[1]
print(V0val)

dum = cc[0].subs(V0,V0val)
print(solve(dum,V1))
