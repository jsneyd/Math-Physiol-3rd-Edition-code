
#     -------------------------------------------------------------------
# 
#      Solve the Topp model for the pathway to diabetes.
# 
#      For Chapter 16, Section 16.7.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Define the ODEs
def rhs(t, y):
    G, I, beta = y
    return [R0 - (EG0 + SI * I) * G,
            beta * sigma * G * G / (alpha * alpha + G * G) - k * I,
            (-d0 + r1 * G - r2 * G * G) * beta]


# Define parameters
SI = 1
EG0 = 1.44
R0 = 864
sigma = 43.2
alpha = 141
k = 432
d0 = 0.06
r1 = 1e-3
r2 = 0.24e-5

# Calculate G1 and G2
G1 = (r1 + np.sqrt(r1 ** 2 - 4 * r2 * d0)) / (2 * r2)
G2 = (r1 - np.sqrt(r1 ** 2 - 4 * r2 * d0)) / (2 * r2)
mu1 = R0 / (alpha * EG0)
Ih1 = (R0 / G1 - EG0) / SI
Ih2 = (R0 / G2 - EG0) / SI
beta1 = (alpha ** 2 + G1 ** 2) / (sigma * G1 ** 2) * k / SI * (R0 / G1 - EG0)
betah1 = k * Ih1 * (alpha ** 2 + G1 ** 2) / G1 ** 2
betah2 = k * Ih2 * (alpha ** 2 + G2 ** 2) / G2 ** 2
mu2 = sigma * SI * betah1 / (EG0 ** 2)
mu3 = k / EG0
eps = d0 / EG0
lam1 = r1 * alpha / d0
lam2 = r2 * alpha ** 2 / d0


# Generate values for G and beta
G = np.linspace(50, 1000, 1000)
beta = (alpha ** 2 + G ** 2) / (sigma * G ** 2) * k / SI * (R0 / G - EG0)

# Plot the phase plane
plt.figure(1)
plt.plot(beta, G, '--', label=r'$\beta$ vs. $G$')
plt.plot([beta1, beta1], [min(G), max(G)], 'k--', label=r'$\beta_1$', linewidth=2)
plt.plot([0, R0 / EG0], [G1, G1], 'k--', label=r'$G_1$', linewidth=2)
plt.plot([0, R0 / EG0], [G2, G2], 'k--', label=r'$G_2$', linewidth=2)
plt.plot([0, 0], [min(G), max(G)], 'k--')
plt.plot([R0 / EG0, R0 / EG0], [min(G), max(G)], 'k--', label=r'$R_0/\varepsilon_0$')
plt.xlim([-5, 100])
plt.xlabel(r'$\beta$')
plt.ylabel('$G$')


# Define the initial conditions and time span
Idataset = [[100, 100, 10], [600, 100, 12], [100, 100, 17], [600, 100, 17]]
tspan = np.linspace(0, 100, 10000)

# Solve the ODEs for each dataset and plot the results
plt.figure(1)
for dataset in Idataset:
    sol = solve_ivp(rhs, [tspan[0], tspan[-1]], dataset, t_eval=tspan)
    plt.plot(sol.y[2], sol.y[0], linewidth=2)
plt.xlabel('beta')
plt.ylabel('G')
plt.xlim([-5, 100])

plt.figure(2)
for dataset in Idataset:
    sol = solve_ivp(rhs, [tspan[0], tspan[-1]], dataset, t_eval=tspan)
    plt.plot(sol.t, sol.y[0])
plt.xlabel('t')
plt.ylabel('G')

plt.tight_layout()
plt.show()
