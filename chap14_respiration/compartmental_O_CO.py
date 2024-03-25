
#   -------------------------------------------------------------------
# 
#   The compartmental model for the removal of carbon monoxide.
# 
#    For Chapter 14, Section 14.4.5 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


# Oxygen saturation function
def conc_O(O, U):
    return O + 4 * THb * O**4 / (KO**4 + O**4 + beta**4 * U**4)


# CO saturation function
def conc_U(O, U):
    return U + 4 * THb * beta**4 * U**4 / (KO**4 + O**4 + beta**4 * U**4)


# ODE for TO and TU
def Tode(t,x):
    global PCO, PO, Otilde, Utilde
    TO = x[0]
    TU = x[1]
    if t > 20:  # Put into a zero CO and high O2 environment after t = 20
        PCO = 0
        PO = PO_new
    [O,Otilde,U,Utilde] = get_O_Otilde_U_Utilde(TO, TU)
    dTO = D * (sigma * PO - Otilde) - gamma * O
    dTU = D * (sigmaCO * PCO - Utilde)
    return [dTO, dTU]


def get_OU(x,TO,TU):
    O = x[0]
    U = x[1]
    return [conc_O(O,U) - TO,  conc_U(O,U) - TU]

def getOtildeUtilde(x,TO,TU):
    global PCO, PO, O, U, Otilde, Utilde
    Otilde = x[0];
    Utilde = x[1];
    OO = D *(sigma*PO - Otilde) + TO - conc_O(Otilde,Utilde);
    UU = D *(sigmaCO*PCO - Utilde) + TU - conc_U(Otilde,Utilde);
    return [OO,UU]

# Function to get O and Otilde as functions of TO and TU
def get_O_Otilde_U_Utilde(TO, TU):
    global PCO, PO, O, U, Otilde, Utilde
    sol_OU = fsolve(lambda x: get_OU(x,TO,TU), [O, U])
    O, U = sol_OU
    sol_tildes = fsolve(lambda x: getOtildeUtilde(x,TO,TU), [Otilde, Utilde])
    Otilde, Utilde = sol_tildes
    return [O,Otilde,U,Utilde]


# One word of warning. Some of the variables must be made global so that, when changed
# one function, the changed values are used in other functions. Same with the initial
# conditions for the solve.

# You have to be very careful with the starting point for the nonlinear solve or else things
# go weird fast.

# Parameters
sigma = 1.4
sigmaCO = 33
beta = 200
KO = 30 * 1.4
THb =2200
RT = 1.7e-2
D = 32.5
PO = 0.21 * (760 - 47)
gamma = 37
PO_newList =  [150, 1 * (760) - 47, 2.5 * (760) - 47]

# Initial conditions
O = 0
U = 0
Otilde = 0
Utilde = 0

# Time array
tdown = 0
tup = 400
tscale = 249/120
tspan = np.linspace(tdown, tup, 1000)


# Solve the ODE for TO and TU for different initial conditions
for ic in range(3):
    PO_new = PO_newList[ic]
    PO = 150
    PCO = 1
    init = [6800, 0]
    sol = solve_ivp(Tode, [tdown,tup], init, t_eval=tspan)
    TO = sol.y[0,:]
    TU = sol.y[1,:]

    plt.figure(ic + 1)
    plt.plot(tspan*tscale,TO, 'r', label='O_2')
    plt.plot(tspan*tscale, TU, 'g', label='CO')
    plt.xlabel('time (minutes)')
    plt.ylabel('total concentration (\u03BCM)')
    plt.title(f'P_O2 = {PO_new} Atm')
    plt.legend()
    M = max(TU)
    half = np.where(TU > M / 2)[0][-1]
    half_clearance = (tspan[half] - 20) * tscale
    plt.show()
    print(f'Half-clearance time for P_O2 = {PO_new} Atm: {half_clearance:.2f} minutes')
