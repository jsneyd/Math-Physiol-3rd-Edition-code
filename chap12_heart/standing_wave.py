
#   -------------------------------------------------------------------
# 
#    Symbolic calculations for the standing wave in a passive cardiac cable.
# 
#    For Chapter 12, Section 12.4.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import sympy as sp

a1, a2, b, L, mu, x, l, ratc, rc, re, rgbyrc, E = sp.symbols('a1 a2 b L mu x l ratc rc re rgbyrc E')

a1 = mu - 1 / E
a2 = mu - E
l = sp.log(E) / L

# Define the equations
ph = a1 * sp.exp(l * x) + a2 * sp.exp(-l * x)
phi = ratc * ph + b
phe = -(1 - ratc) * ph + b

eq1 = phe.subs(x, L) - mu * phe.subs(x, 0)

# Solve for b
sol = sp.solve(eq1, b)
solb = sp.simplify(sol[0])

phi = phi.subs(b, solb)
dphi = sp.diff(phi, x)
dphiL = dphi.subs(x, L)
eq4 = phi.subs(x, L) - mu * phi.subs(x, 0) + rgbyrc * dphiL

print(sp.pretty(sp.simplify(eq4)))   # This is the equation to solve for mu
