#   -------------------------------------------------------------------
# 
#    Program to compute the response of the AMPA model to a delta function
#    input
# 
#    For Chapter 8, Fig. 8.13 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the system of ODEs
def rhs(t,s):
    c0, c1, c2, D1, D2, simpleC = s
    O = 1 - c0 - c1 - c2 - D1 - D2
    TT = 1000*(np.heaviside(t,0) - np.heaviside(t-0.001,0))

    c0p = -Rb * TT * c0 + Ru1 * c1
    c1p = Rb * TT * c0 - Ru1 * c1 + Rr * D1 - Rd * c1 + Ru2 * c2 - Rb * TT * c1
    c2p = -Ru2 * c2 + Rb * TT * c1 + Rr * D2 - Rd * c2 + Rc * O - Ro * c2
    D1p = Rd * c1 - Rr * D1
    D2p = Rd * c2 - Rr * D2
    simpleCp = -alpha * simpleC * TT + beta * (1 - simpleC)

    return [c0p, c1p, c2p, D1p, D2p, simpleCp]


# Parameters
Rb = 13.0e6
Ru1 = 5.9
Ru2 = 8.6e4
Rd = 900.0
Rr = 64.0
Ro = 2.7e3
Rc = 200.0
alpha = 1.35
beta = 200.0

# Time settings
t_end = 0.2
tout= np.linspace(0,t_end,1000)

# Initial conditions
s0 = [1, 0, 0, 0, 0, 1]

# Solve the system of ODEs. LSODA doesn't work. Have to use BDF.
solution = solve_ivp(rhs, [0,t_end], s0, t_eval=tout, method='BDF')

# Extract results
c0 = solution.y[0]
c1 = solution.y[1]
c2 = solution.y[2]
D1 = solution.y[3]
D2 = solution.y[4]
simpleC = solution.y[5]
simpleO = 1 - simpleC
O = 1 - c0 - c1 - c2 - D1 - D2

# Plot results
plt.plot(solution.t, O, label='Full Model')
plt.plot(solution.t, simpleO, '--', label='Simple Model')
plt.xlabel('time (ms)')
plt.ylabel('Fraction of Open Channels')
plt.legend()
plt.show()
