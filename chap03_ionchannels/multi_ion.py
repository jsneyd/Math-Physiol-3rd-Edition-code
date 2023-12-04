# -------------------------------------------------------------------

#  Python code to find analytic expressions for the flux through a
#  saturating barrier model of an ion channel, with N binding sites. Here,
#  the channel can bind two ions at once.

#  For Chapter 3, Section 1.4.3 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

from sympy import symbols, solve, fraction, simplify, limit, pretty

k12, k13, k14, k21, k23, k24, k31, k32, k34, k41, k42, k43, kf, kon = symbols('k12 k13 k14 k21 k23 k24 k31 k32 k34 k41 k42 k43 kf kon')
p1, p2, p3, p4, ci, ce = symbols('p1 p2 p3 p4 ci ce')

a1 = 2*kf
a2 = kon*ce + kf + k42
a3 = kon*ci + kon*ce
a4 = kon*ci + k24 + kf

eq1 = -a1*p1 + kon*ce*p2 + kon*ci*p4
eq2 = kf*p1 - a2*p2 + kon*ci*p3 + k24*p4
eq3 = kf*p2 - a3*p3 + kf*p4
eq4 = p1 + p2 + p3 + p4 - 1

sols = solve([eq1, eq2, eq3, eq4], [p1, p2, p3, p4])
p2 = sols[p2]
p4 = sols[p4]

J = p2*k42 - p4*k24

# Display the numerator of J
NJ, DJ = J.as_numer_denom()
print(pretty(simplify(NJ)))

# Calculate limits
limit_kf_inf = limit(J, kf, float('inf'))
limit_ci_inf = limit(J, ci, float('inf'))
limit_ce_inf = limit(J, ce, float('inf'))

print("Limit as kf approaches infinity:", limit_kf_inf)
print("Limit as ci approaches infinity:", limit_ci_inf)
print("Limit as ce approaches infinity:", limit_ce_inf)