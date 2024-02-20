
#   -------------------------------------------------------------------
# 
#    Symbolic computations for a three-compartment model of the circulatory
#    system.
# 
#    For Chapter 11, Section 11.5.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

from sympy import symbols, solve, simplify

# Define symbolic variables
Pa, Pv, Va, Vv, Q = symbols('Pa Pv Va Vv Q')
Cd, F, Rs, Ca, Cv, Vt = symbols('Cd F Rs Ca Cv Vt')

# Define symbolic equations
fun1 = Q - F*Cd*Pv
fun2 = Q - (Pa - Pv)/Rs
fun3 = Va - Ca*Pa
fun4 = Vv - Cv*Pv

# Volume equation
fun5 = Vt - Va - Vv

# Define the system of equations
eqns = [fun1, fun2, fun3, fun4, fun5]

# Solve the equations
sols = solve(eqns, [Pa, Pv, Va, Vv, Q])

# Simplify the solution for Q
Q = simplify(sols[Q], rational=True)
print(Q)