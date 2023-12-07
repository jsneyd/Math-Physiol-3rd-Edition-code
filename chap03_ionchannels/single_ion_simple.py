# -------------------------------------------------------------------

#  Python code to find analytic expressions for the flux through a
#  saturating barrier model of an ion channel, with N binding sites, where
#  every barrier is symmetric and identical.

#  For Chapter 3, Exercise 3.6 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import sympy as sp

x, ci, ce, kon, J, k, km, konL, konR, v, kbar = \
        sp.symbols('x ci ce kon J k km konL konR v kbar')

# N is the number of binding sites.
N = 2

eq = [sp.symbols('eq{}'.format(i)) for i in range(1, N)]
c = [sp.symbols('c{}'.format(i)) for i in range(1, N+1)]

# Define the endpoint equations
eq0 = sp.Eq(konL*ci*x - km*c[0], J)
eqNN = sp.Eq(k*c[-1] - konR*ce*x, J)

# Define the internal equations
for i in range(N-1):
    eq[i] = sp.Eq(k*c[i] - km*c[i+1], J)

# Define the conservation equation
eqcons = sp.Eq(x + sum(c), 1)


# Solve for J in terms of c
sol = sp.solve([eq0] + eq + [eqNN] + [eqcons], [J]+c+[x])
print(sp.pretty(sp.simplify(sol[J])))


# Make the rate constants functions of v = VF/RT

subs_V = [(konL, kon * sp.exp(v / (2 * (N + 1)))),
         (konR, kon * sp.exp(-v / (2 * (N + 1)))),
         (k, kbar * sp.exp(v / (2 * (N + 1))) ),
         (km, kbar * sp.exp(-v / (2 * (N + 1))))]        
             
# Substitute the expressions for k0, k, and km into the solution for J
J_subs = sol[J].subs(subs_V)
print(sp.pretty(sp.simplify(J_subs)))


print('\n\n max flux', (sp.simplify(sp.limit(J_subs,ci,sp.oo))),'\n\n')


## This next bit is for calculating the derivative of Jmax

down = 1

# Calculate the sum in the loop
for i in range(2, N + 1):
    down = down + i * sp.exp((i - 1) * v / (N + 1))

# Calculate J_max
J_max = sp.exp((2 * N - 1) * v / (2 * N + 2)) / down

# Calculate and display the derivative of J_max with respect to v
derivative_J_max = sp.diff(J_max, v)
print( sp.pretty(sp.simplify(derivative_J_max)))


