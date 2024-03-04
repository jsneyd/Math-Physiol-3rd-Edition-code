import sympy as sp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


# First do the symbolic calcuation of the steady state, and plot it

k1, km1, k2, km2, R, O, I, c = sp.symbols('k1 km1 k2 km2 R O I c')

RI = 1 - R - O - I
eq1 = km1*O + km2*RI - R*(k1*c*c + k2*c)
eq2 = km2*I + k1*c*c*R - O*(k2*c + km1)
eq3 = k2*c*O + k1*c*c*RI - I*(km2 + km1)

k1val = 35
km1val = 0.06
k2val = 0.5
km2val = 0.01

sols = sp.solve([eq1, eq2, eq3], (R, O, I))
Ost = sp.lambdify(c, sols[O].subs([(k1, k1val), (km1, km1val), (k2, k2val), (km2, km2val)]))
Rst = sp.lambdify(c, sols[R].subs([(k1, k1val), (km1, km1val), (k2, k2val), (km2, km2val)]))
Ist = sp.lambdify(c, sols[I].subs([(k1, k1val), (km1, km1val), (k2, k2val), (km2, km2val)]))

plt.figure(1)
c_vals = np.linspace(0.01, 0.5, 100)
plt.plot(c_vals, Rst(c_vals), 'r', label='R')
plt.plot(c_vals, Ost(c_vals), 'b', label='O (steady)')
plt.plot(c_vals, Ist(c_vals), 'g', label='I')
plt.xlabel('c')




# Next, do the numerical solutions and calculate the peak response

# Define the ODEs
def rhs(t, x, c):
    global k1val, km1val, k2val, km2val
    R, O, I = x
    RI = 1 - R - O - I

    k1 = k1val
    km1 = km1val
    k2 = k2val
    km2 = km2val

    return [
        km1 * O + km2 * RI - R * (k1 * c * c + k2 * c),
        km2 * I + k1 * c * c * R - O * (k2 * c + km1),
        k2 * c * O + k1 * c * c * RI - I * (km2 + km1)
    ]

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

# Calculate the steady-state and peak responses
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


