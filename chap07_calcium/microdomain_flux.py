#  code to solve the microdomain flux equations

#    For Chapter 7, Section 7.6 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt



# Function to compute w_v
def getw_v(x, D):
    return D * x + Dbbt * x / (Kb + x)

# Function to compute the derivative of w_v
def getw_v_prime(x, D):
    return D + Dbbt * Kb / (Kb + x)**2

# Function to compute ce0
def getce0(c0):
    w0 = getw_v(c0, Dc)
    wp0 = getw_v_prime(c0, Dc)
    return c0 + Dc * (w0 - winf) / (wp0 * rho)

# Function for c0 that must be set to zero to find c0
def getc0(c0):
    ce0 = getce0(c0)
    w0 = getw_v(c0, Dc)
    wp0 = getw_v_prime(c0, Dc)
    v0 = getw_v(ce0, De)
    vp0 = getw_v_prime(ce0, De)
    return Dc * (w0 - winf) / wp0 + De * (v0 - vinf) / vp0


# Global variables
Dc = 1
De = 1
Dbbt = 5
rho = 2
cinf = 0.5

ceinf_list = [1, 5, 10]
n = 1500
Kb_list = np.linspace(0.001, 10, n)
rho_eff = np.zeros(n)  # preallocate space

# Loop over ceinf_list and Kb_list to calculate rho_eff and plot
for ceinf in ceinf_list:
    for i in range(n):
        Kb = Kb_list[i]
        winf = getw_v(cinf, Dc)
        vinf = getw_v(ceinf, De)
        c0_initial_guess = 2  # initial guess for fsolve
        c0 = fsolve(getc0, c0_initial_guess)
        J = rho * (getce0(c0) - c0)
        rho_eff[i] = J / (ceinf - cinf)
    
    plt.plot(Kb_list, rho_eff, label=f'c_e,∞ = {ceinf}')

# Add labels and legend to the plot
plt.xlabel('K_b')
plt.ylabel('ρ_eff')
plt.legend()
plt.show()
