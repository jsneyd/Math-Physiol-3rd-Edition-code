#   -------------------------------------------------------------------
#    Shooting to find spiral wave profile
#
#    For Chapter 6, Figure 6.15 and 6.16,   of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


def deRHS(r, s):
    ph = s[0]
    th = s[1]
    ps = ph / np.sqrt(r**2 - ph**2)
    dph = r * (c - w * np.sqrt(r**2 - ph**2))
    dth = ps / r
    return [dph, dth]

def events(t, s):
    value = abs(s) - 100
    isterminal = 1
    direction = 0
    return value, isterminal, direction

# Set plot parameters
plt.rcParams.update({
    'axes.titlesize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0,
    'patch.linewidth': 0.7
})



# Global variables
c = 3.0
w = 0
eps = 0.05
alp = 0.1
mstar = 0.3318
wlist = [2.8, 9 * mstar, 3.2]
for j in range(3):
    w = wlist[j]
    rstep = -0.001
    r_end = 10
    s0 = [r_end - 0.1, 0]
    tspan = np.arange(r_end, 0, rstep)
    
    sol = solve_ivp(deRHS, [r_end, 0], s0, t_eval=tspan, method='RK45')
    
    r = sol.t
    S = sol.y.T
    ph = S[:, 0]
    th = S[:, 1]

    ps = ph / np.sqrt(r**2 - ph**2)
    plt.figure(1)
    ndx = np.where(ph**2 < r**2)[0]
    plt.plot(r[ndx], ps[ndx])
    plt.xlabel('r')
    plt.ylabel(r'$\psi$')
    plt.axis([0, 5, -4.5, 4.5])
    
    if j == 2:
        r0 = c / w
        R = np.arange(r0, r_end, 0.01)
        rho = np.sqrt(R**2 / r0**2 - 1)
        th0 = 2 * np.pi * np.arange(0, 1.01, 0.01)
        X0 = r0 * np.cos(th0)
        Y0 = r0 * np.sin(th0)
        s = rho - rho[-1] + th[0] + np.pi / 2
        X1 = r0 * np.cos(s) + r0 * rho * np.sin(s)
        Y1 = r0 * np.sin(s) - r0 * rho * np.cos(s)

        plt.figure(2)
        X = r * np.cos(th)
        Y = r * np.sin(th)
        plt.plot(X, Y, X0, Y0, '--', X1, Y1, '--')
        plt.axis([-10, 10, -10, 10])

plt.figure(1)
plt.legend(['ω/ε=2.8', 'ω/ε=2.9864', 'ω/ε=3.2'])

# Singular dispersion curve
w0 = np.arange(0, (1 - 2 * alp) / 2, 0.001)
apw = alp + w0
w1 = 1 - 2 * alp - w0
sp = (1 - 2 * apw) / np.sqrt(apw - apw**2)
T = np.log(((1 - w0) * w1) / ((1 - w1) * w0))
ccrit = np.sqrt(eps / (T * mstar))

plt.figure(10)
plt.plot(2 * np.pi / T, sp, '--', 1 / T, ccrit)
plt.axis([0, 10, 0, 3])
plt.xlabel('Frequency, 2π/T')
plt.ylabel('Speed, c')
plt.text(1.8, 2.3, 'dispersion curve', fontsize=18)
plt.text(7, 0.9, 'critical curve', fontsize=18)


