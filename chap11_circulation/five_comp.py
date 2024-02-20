
#   -------------------------------------------------------------------
# 
#    Symbolic computations for a five-compartment model of the circulatory
#    system.
# 
#    For Chapter 11, Section 11.5.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

from sympy import symbols, solve, simplify, limit, pretty

# Define symbolic variables
QQ, Psa, Psv, Ps, Ppa, Ppv, Vsa, Vsv, Vp = symbols('QQ Psa Psv Ps Ppa Ppv Vsa Vsv Vp')
F, Cld, Crd, Csa, Csv, Cp, Rsa, Rsv, Rp, Vt, V0p, V0s = symbols('F Cld Crd Csa Csv Cp Rsa Rsv Rp Vt V0p V0s')

# Systemic arteries
fun0 = QQ - (Psa - Ps)/Rsa
fun1 = QQ - F*Cld*Ppv
fun2 = Vsa - V0s - Csa*(Psa + Ps)/2

# Systemic veins
fun3 = QQ - (Ps - Psv)/Rsv
fun4 = Vsv - Csv*(Psv + Ps)/2

# Pulmonary system
fun4a = QQ - (Ppa - Ppv)/Rp
fun5 = QQ - F*Crd*Psv
fun6 = Vp - V0p - Cp*(Ppa + Ppv)/2

# Volume
fun7 = Vt - Vsa - Vsv - Vp

# Define the system of equations
eqns = [fun0, fun1, fun2, fun3, fun4, fun4a, fun5, fun6, fun7]

# Solve the equations
sols = solve(eqns, [Psa, Psv, Ps, Ppa, Ppv, Vsa, Vsv, Vp, QQ])

# Simplify the solution for QQ
QQ = simplify(sols[QQ], rational=True)

# Compute the limit as F approaches infinity
Qinf = simplify(limit(QQ, F, float('inf')))

# Print the result in a pretty format
print(pretty(Qinf))
