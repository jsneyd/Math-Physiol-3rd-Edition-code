
#   -------------------------------------------------------------------
# 
#    Phase portraits for leukocyte dynamics. With non-zero alpha.
# 
#    For Chapter 13, Section 13.2.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def rhs(x, t):
    global N
    out = np.zeros(2 * N)
    for j in range(N):
        u = x[2*j]
        v = x[2*j + 1]
        f = F(a * v)
        up = g0 * ((b0 + b * v) * (1 - u * f * np.exp(-a * v)) - u * (v * f + 1))
        vp = v * (1 - xi * u * f)
        out[2*j] = up
        out[2*j + 1] = vp
    return out


def F(z):
    z = np.array(z)             # don't ask why, you just have to
    z[z<0.00001] = 0.00001      # avoid singularity at z=0
    f = z / (1 - np.exp(-z))
    return f


# Parameters
g0 = 0.2
b0 = 1
a = 0

xivals = [1.6, 3.0, 1.6, 3]
bvals = [0.1, 0.1, 3, 3]
letters = ['A', 'B', 'C', 'D']
vinit = [0.5, 2.0, 1.5, 1.5, 2.1]
U0 = b0 / (1 + b0)

for j in range(4):
    b = bvals[j]
    xi = xivals[j]
    letter = letters[j]
    v = np.linspace(0, 2.5, 10) #227)
    f = F(a * v)
    u1 = (b0 + b * v) / ((b0 + b * v) * f * np.exp(-a * v) + (v * f + 1))
    u2 = 1 / (xi * f)
    diff = u1 - u2
    ndx = np.where(diff[:-1] * diff[1:] <= 0)[0]

    N = 1
    init = [U0, vinit[j]]
    if j == 1:
        N = 2
        init = [U0, vinit[j], U0, vinit[4]]

    dt = 0.01
    tend = [5, 30, 30, 30]
    tspan = np.arange(0, tend[j] + dt, dt)
    sol = odeint(rhs, init, tspan)

    plt.figure(j+1)
    plt.plot(sol[:, 0], sol[:, 1], 'b', linewidth=2)
    plt.plot(u1, v, 'g--', u2, v, 'r--', linewidth=2)
    plt.plot(U0, 0, '*', linewidth=2)
    if len(ndx) > 0:
        plt.plot(u2[ndx], v[ndx], '*', linewidth=2)
    plt.xlabel('U')
    plt.ylabel('V')
    plt.title('{}: $\\xi = {:.2f}$, $\\beta = {:.2f}$'.format(letters[j], xi, b))
    plt.axis([0.3, 0.7, 0, 2.5])
    plt.grid(True)


# Plot parameter curves
V = np.arange(0.1, 1000.1, 0.1)
alist = [0.04, 0.5, 0.75, 1]

for a in alist:
    cA = a**2 * V**4
    cB = -V * (-2 * a**2 * V**2 * b0 + a * V - 2 * np.exp(a * V) + 2)
    cC = b0 * (a**2 * V**2 * b0 - a * V - 1 + np.exp(a * V) * (1 - V**2 * a))
    disc = cB**2 - 4 * cA * cC
    b1 = (-cB + np.sqrt(disc)) / (2 * cA)

    ndx = np.where(b1 > 0)
    onebyxi1 = a * V * (b0 + b1 * V) / (1 - np.exp(-a * V) + a * V**2 + a * V * np.exp(-a * V) * (b0 + b1 * V))

    plt.figure(6)
    plt.plot(1 / b1[ndx], 1 / onebyxi1[ndx],label=r'$\alpha = ${}'.format(a))

plt.xlabel('1/β')
plt.ylabel('ξ')


plt.figure(6)
plt.plot([0, 4], [0, 4], '--', [0, 4], [1 / U0, 1 / U0], '--')
plt.xlabel(r'$1/\beta$')
plt.ylabel(r'$\xi$')
plt.axis([0,4,0,4])
plt.legend()





