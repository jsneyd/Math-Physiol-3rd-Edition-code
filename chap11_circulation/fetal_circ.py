
#   -------------------------------------------------------------------
# 
#    The model of fetal circulation.
# 
#    For Chapter 11, Section 11.7 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import sympy as sp

# Define the symbolic variables
Psa, Psv, Ps, Rsa, Rsv, Ppa, Ppv = sp.symbols('Psa Psv Ps Rsa Rsv Ppa Ppv')
Cld, Crd, Csv, Csa, Cp = sp.symbols('Cld Crd Csv Csa Cp')
Vsv, V0p, Vp, V0s, Vsa, Vt = sp.symbols('Vsv V0p Vp V0s Vsa Vt')
Rd, Rp, F = sp.symbols('Rd Rp F')
Qr, Qs, Qd, Qp, Ql, Qf = sp.symbols('Qr Qs Qd Qp Ql Qf')

# just to make things simpler
Crs = 0
Cls = 0

# initialise
eq = sp.zeros(13)
eq_open = sp.zeros(13)
eq_closed = sp.zeros(13)

# define the equations
eq[0] = Psa - Ps - Qs * Rsa
eq[1] = Ps - Psv - Qs * Rsv
eq[2] = Ppa - Ppv - Qp * Rp
eq[3] = F * (Cld * Ppv - Cls * Psa) - Ql
eq[4] = F * (Crd * Psv - Crs * Ppa) - Qr
eq[5] = Ql + Qd - Qs
eq[6] = Qd + Qp - Qr
eq[7] = Qr + Qf - Qs
eq[8] = Csv * (Psv + Ps) / 2 - Vsv
eq[9] = Cp * (Ppa + Ppv) / 2 + V0p - Vp
eq[10] = Csa * (Psa + Ps) / 2 + V0s - Vsa
eq[11] = Ppa - Psa - Qd * Rd
eq[12] = Vt - Vsa - Vsv - Vp


# The foramen-open solution
print('The foramen-open solution\n\n')
# set Psv = Ppv
for j in range(13):
    eq_open[j] = eq[j].subs(Psv, Ppv)

sol_open = sp.solve(eq_open, 
                    [Psa, Ppa, Ppv, Ps, Vp, Vsv, Vsa, Vt, Qf, Qs, Qd, Qp, Ql])

# Print selected solution
print(sp.pretty(sp.simplify(sol_open[Qs])),'\n')
print(sp.pretty(sp.simplify(sol_open[Qp])),'\n')
print(sp.pretty(sp.simplify(sol_open[Ql])),'\n')


# The foramen-closed solution
print('\n\nThe foramen-open solution\n\n')
# set Qf = 0
for j in range(13):
    eq_closed[j] = eq[j].subs(Qf, 0)

sol_closed = sp.solve(eq_closed, 
                    [Psa, Ppa, Ppv, Psv, Ps, Vp, Vsv, Vsa, Vt, Qs, Qd, Qp, Ql])

# Print selected solution
print(sp.pretty(sp.simplify(sol_closed[Ql])),'\n')
print(sp.pretty(sp.simplify(sol_closed[Ppv])),'\n')
print(sp.pretty(sp.simplify( sol_closed[Ppv]/sol_closed[Psv] )))
