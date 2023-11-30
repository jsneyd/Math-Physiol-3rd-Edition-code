# -------------------------------------------------------------------

#  Python code to find analytic expressions for the flux through a
#  saturating barrier model of an ion channel, with N binding sites

#  For Chapter 3, Section 1.4.2 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

from sympy import symbols, Eq, solve, pretty, simplify

x, ci, ce, k0, J, kNN, kmNN = symbols('x ci ce k0 J kNN kmNN')

# N is the number of binding sites.
N = 1
eq = [symbols('eq{}'.format(i)) for i in range(1, N)]

k = [symbols('k{}'.format(i)) for i in range(1, N+1)]
km = [symbols('km{}'.format(i)) for i in range(1,N+2)]
c = [symbols('c{}'.format(i)) for i in range(1, N+1)]

# Define the endpoint equations
eq0 = Eq(k0*ci*x - km[0]*c[0], J)
eqNN = Eq(k[-1]*c[-1] - km[-1]*ce*x, J)

# Define the internal equations
for i in range(N-1):
    eq[i] = Eq(k[i]*c[i] - km[i+1]*c[i+1], J)

# Define the conservation equation
eqcons = Eq(x + sum(c), 1)

# Solve for J in terms of x
sol_J = solve([eq0] + eq + [eqNN], [J]+c)
pretty_sol_J = pretty(simplify(sol_J[J]))
print(pretty_sol_J)

# Solve for J in terms of x, c
sol = solve([eq0] + eq + [eqNN] + [eqcons], c + [J])
pretty_sol_J_cons = pretty(simplify(sol[J]))
print(pretty_sol_J_cons)
