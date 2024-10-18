#   -------------------------------------------------------------------
# 
#    Reponse of the RyR model of Stern to a step increase in calcium.
# 
#    For Chapter 7, Section 7.4.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------


import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


# Calculate the steady states for any given c
def ROIss(c):
    A = np.array([[-km2 - (k1 * c**2 + k2 * c), km1 - km2, -km2],
                  [k1 * c**2, -(k2 * c + km1), km2],
                  [-k1 * c**2, -k1 * c**2 + k2 * c, -(km2 + km1) - k1 * c**2]])
    
    Rhs = np.array([-km2, 0, -k1 * c**2])
    R = np.linalg.solve(A, Rhs)
    return R[0], R[1], R[2]


# Define the ODEs
def rhs(t, x, c):
    R, O, I = x
    RI = 1 - R - O - I

    return [
        km1 * O + km2 * RI - R * (k1 * c * c + k2 * c),
        km2 * I + k1 * c * c * R - O * (k2 * c + km1),
        k2 * c * O + k1 * c * c * RI - I * (km2 + km1)
    ]


k1 = 35
km1 = 0.06
k2 = 0.5
km2 = 0.01

# Plot the steady states as functions of c
clist = np.linspace(0.01, 0.5, 50)
Rst, Ost, Ist = [], [], []
for c in clist:
    R, O, I = ROIss(c)
    Rst.append(R)
    Ost.append(O)
    Ist.append(I)

plt.figure(1)
plt.plot(clist, Rst, 'r', label='Rst')
plt.plot(clist, Ost, 'b', label='Ost')
plt.plot(clist, Ist, 'g', label='Ist')
plt.legend()
plt.xlabel('c (\u03BCM)')
plt.ylabel('Steady states')
plt.title('Steady states vs. c')



# Next, do the numerical solutions and calculate the peak response

# Set initial conditions, time span, and parameter values
init = [1, 0, 0]  # Initial conditions for R, O, I
tspan = np.linspace(0, 1.5, 1000)  # Time span
c_vals = [1, 2, 3]  # Values of c to solve for

# Solve the ODEs for different values of c
plt.figure(2)
for c in c_vals:
    sol = solve_ivp(lambda t, x: rhs(t, x, c), (tspan[0], tspan[-1]), init, t_eval=tspan)
    plt.plot(sol.t, sol.y[1], label=f'c = {c}')

plt.xlabel('Time')
plt.ylabel('O')
plt.legend()

# Calculate and plot the peak response
cc = np.linspace(0.01, 0.5, 50)
tspan = np.linspace(0, 50, 100000)
peak = []
for c in cc:
    sol = solve_ivp(lambda t, x: rhs(t, x, c), (tspan[0], tspan[-1]), init, t_eval=tspan)
    peak.append(sol.y[1].max())

plt.figure(1)
plt.plot(cc, peak, 'k',label='O (peak)')
plt.xlabel('c')
plt.legend()


