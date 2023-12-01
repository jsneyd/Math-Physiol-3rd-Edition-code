# -------------------------------------------------------------------

#  Python code to find analytic expressions for the flux through a
#  nonsaturating barrier model of an ion channel, with N binding sites

#  For Chapter 3, Section 1.4.1 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

from sympy import symbols, solve, pretty, simplify

x, ci, ce, k0, M = symbols('x ci ce k0 M')

# N is the number of binding sites. Note that with N binding sites you have 2
# endpoint equations and an additional N-1 equations for the binding sites.

# The text is written in terms of N barriers, and thus N-1 binding sites,
# but it's slightly easier to write the code for N binding sites. 


N = 3
eq = [symbols('eq{}'.format(i)) for i in range(1, N)]

k = [symbols('k{}'.format(i)) for i in range(1, N+1)]
km = [symbols('km{}'.format(i)) for i in range(1,N+2)]
c = [symbols('c{}'.format(i)) for i in range(1, N+1)]

# Define the endpoint equations
eq0 = k0*ci - km[0]*c[0] - M
eqNN = k[-1]*c[-1] - km[-1]*ce -  M

# Define the internal equations
for i in range(N-1):
    eq[i] = k[i]*c[i] - km[i+1]*c[i+1] - M
    print(i)

# Solve for J
sol = solve([eq0]+eq+[eqNN],c+[M])
result_M = sol[M]

# Display the result
pretty_result_M = pretty(simplify(result_M))
print(pretty_result_M)
