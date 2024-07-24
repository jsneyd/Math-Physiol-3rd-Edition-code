# -------------------------------------------------------------------

#  Python code to find analytic expressions for the flux through a
#  saturating barrier model of an ion channel, with 2 binding sites. Here,
#  the channel can bind two ions at once.

#  For Chapter 3, Section 3.4.3 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import sympy as sp

# Define symbols
k0, k1, k2, km1, km2, km3 = sp.symbols('k0 k1 k2 km1 km2 km3')
p00, p01, p10, p11 = sp.symbols('p00 p01 p10 p11')
ci, ce = sp.symbols('ci ce')

# Define equations
a1 = k0 * ci + km3 * ce
a2 = km1 + k1 + km3 * ce
a3 = k2 + km2 + k0 * ci
a4 = k2 + km1

eq1 = sp.Eq(-a1 * p00 + km1 * p10 + k2 * p01, 0)
eq2 = sp.Eq(k0 * ci * p00 - a2 * p10 + km2 * p01 + k2 * p11, 0)
eq3 = sp.Eq(km3 * ce * p00 + k1 * p10 - a3 * p01 + km1 * p11, 0)
eq4 = sp.Eq(p00 + p01 + p10 + p11, 1)

# Solve the system of equations
sols = sp.solve([eq1, eq2, eq3, eq4], (p00, p01, p10, p11))
p10_sol = sols[p10]
p01_sol = sols[p01]

# Define the flux J
J = p10_sol * k1 - p01_sol * km2

# Get the numerator and denominator of J
NJ, DJ = J.as_numer_denom()

# Print the simplified numerator and denominator
print("Numerator:")
sp.pretty_print(sp.factor(NJ))
print("\nDenominator:")
sp.pretty_print(sp.simplify(DJ))

# Limits of J
J_ci_inf = sp.simplify(sp.limit(J, ci, sp.oo))
J_ce_inf = sp.simplify(sp.limit(J, ce, sp.oo))

print("\nLimit of J as ci -> ∞:")
sp.pretty_print(J_ci_inf)
print("\nLimit of J as ce -> ∞:")
sp.pretty_print(J_ce_inf)