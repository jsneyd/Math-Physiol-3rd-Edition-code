# -------------------------------------------------------------------

#  Python code to find analytic expressions for the flux through a
#  saturating barrier model of an ion channel, with N binding sites

#  For Chapter 3, Section 1.4.2 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import sympy as sp

x, ci, ce, k0, kap0, J, kNN, kmNN, v = sp.symbols('x ci ce k0 kap0 J kNN kmNN v')

# N is the number of binding sites.
N = 1

eq = [sp.symbols('eq{}'.format(i)) for i in range(1, N)]
k = [sp.symbols('k{}'.format(i)) for i in range(1, N+1)]
kap = [sp.symbols('kap{}'.format(i)) for i in range(1, N+1)]
km = [sp.symbols('km{}'.format(i)) for i in range(1,N+2)]
kapm = [sp.symbols('kapm{}'.format(i)) for i in range(1,N+2)]
c = [sp.symbols('c{}'.format(i)) for i in range(1, N+1)]

# Define the endpoint equations
eq0 = sp.Eq(k0*ci*x - km[0]*c[0], J)
eqNN = sp.Eq(k[-1]*c[-1] - km[-1]*ce*x, J)

# Define the internal equations
for i in range(N-1):
    eq[i] = sp.Eq(k[i]*c[i] - km[i+1]*c[i+1], J)

# Define the conservation equation
eqcons = sp.Eq(x + sum(c), 1)

# Solve for J in terms of x and c
sol_J = sp.solve([eq0] + eq + [eqNN], [J]+c)
print(sp.pretty(sp.simplify(sol_J[J])))


# Solve for J in terms of c
sol = sp.solve([eq0] + eq + [eqNN] + [eqcons], [J]+c+[x])
print(sp.pretty(sp.simplify(sol[J])))


# Make the rate constants functions of v = VF/RT
# This is a weird, weird construction due to John Rugis. A perfect example 
# of how Python, although good at a lot, really falls down on some
# things.

# If you can think of a better way to do this substitution, do please tell us.

subs_V = [(k0, kap0 * sp.exp(v / (2 * (N + 1))))] + \
          list(zip(k, [j * sp.exp(v / (2 * (N + 1))) for j in kap])) + \
          list(zip(km, [j* sp.exp(v / (-2 * (N + 1))) for j in kapm])) 
             
             
# Substitute the expressions for k0, k, and km into the solution for J
J_subs = sol[J].subs(subs_V)
print(sp.pretty(sp.simplify(J_subs)))
