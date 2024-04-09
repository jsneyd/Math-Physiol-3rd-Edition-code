
#     -------------------------------------------------------------------
# 
#  Solve the Lacker model for the control of ovulation number. The model is
#  solved in the transformed coordinates (gamma rather than xi) because it is
#  much easier that way, as this avoids the problem of  solutions going to infinity in
#  finite time, which is a real pain
# 
#      For Chapter 16, Section 16.3.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


# Define the right-hand side of the ODE system
def rhs(t, y):
    gam = y[:n]
    xi1 = y[n]
    GG = np.sum(gam)
    gam_deriv = gam * (1 - gam) * (M1 * M2 * (1 + gam) - GG * (M1 + M2))
    xi_deriv = (1 / xi1) * (1 - xi1 ** 2 * (GG - M1) * (GG - M2))
    time_deriv = 1 / (xi1 ** 2)
    return np.concatenate((gam_deriv, [xi_deriv, time_deriv]))


# Define parameters
n = 10
M1 = 4
M2 = 20

# Initial conditions
fol0 = np.sort(0.1 * np.random.rand(n))[::-1]  # They must be ordered.
xi10 = fol0[0]
gam0 = fol0 / fol0[0]  # Initial conditions for the gammas.
y0 = np.concatenate((gam0, [xi10, 0]))

# Define the time span
tspan = [0, 2]
tout=np.linspace(tspan[0], tspan[1], 1000)
# Solve the ODE system
sol = solve_ivp(rhs, tspan, y0, method='LSODA', t_eval=tout)

# Extract the results
gam = sol.y[:n]
xi1 = sol.y[n]
t = sol.y[n + 1]

# Plot the results
plt.figure(1)
plt.plot(t, gam.T)
plt.xlabel('tau')
plt.ylabel('gamma_i')

ov_index = np.argmax(xi1 > 5)
t_ov = t[ov_index]
ov_num = np.sum(gam[:, -1] > 0.9)
print("Ovulation time:", t_ov)
print("Number of ovulations:", ov_num)

plt.figure(2)
plt.plot(t, xi1 * gam.T)
plt.xlabel('t')
plt.ylabel('xi_i')
plt.ylim([0, 5])

plt.show()
