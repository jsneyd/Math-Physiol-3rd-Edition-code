
#   -------------------------------------------------------------------
# 
#    Autoregulation in the three-compartment circulation model
# 
#    For Chapter 11, Section 11.6.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import sympy as sp

# Define symbolic variables
Pa, Pv, Va, Vv, Q, Oa, Ov = sp.symbols('Pa Pv Va Vv Q Oa Ov')
Cd, F, Ca, Cv, Vt, Rs = sp.symbols('Cd F Ca Cv Vt Rs')
Smax, beta, Pastar, M, A, KP = sp.symbols('Smax beta Pastar M A KP')

# Define autoregulation function
S = KP * Smax / (KP + Pa)
Rs = S + A * Ov

# Define equations
fun1 = Q - F * Cd * Pv
fun2 = Q - (Pa - Pv) / Rs
fun3 = Va - Ca * Pa
fun4 = Vv - Cv * Pv
fun5 = Vt - Va - Vv             # volume
afun1 = M - Q * (Oa - Ov)       # autoregulation

# Solve equations
sols = sp.solve((fun1, fun2, fun3, fun4, fun5, afun1), (F, Pv, Va, Vv, Q, Ov))

Q = sp.simplify(sols[0][4])
print(sp.pretty(sp.simplify(Q)))

F = sp.simplify(sols[0][0])
print(sp.pretty(sp.simplify(F)))


# remove Vt by substitution. Just because we can. It's not important.
Q = Q.subs(Vt,Pa*Ca+Pv*Cv)
print(sp.pretty(sp.simplify(Q)))

F = F.subs(Vt,Pa*Ca+Pv*Cv)
print(sp.pretty(sp.simplify(F)))


