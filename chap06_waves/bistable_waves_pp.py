#   -------------------------------------------------------------------
# 
#   Plot some solutions of the bistable equation for different values of 
#   c. A demonstration of how shooting can work to find the wave speed.
#
#    The wave speed is slightly different from that calculated by the Matlab code, 
#    but closer to the correct value of (1-2*alp)/(sqrt(2)) = 0.56568542.
# 
#    For Chapter 6, Figure 6.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd
# 
#   -------------------------------------------------------------------

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def de_rhs(t, s):
    v = s[0]
    w = s[1]

    Fv = w
    Fw = c * w - v * (1 - v) * (v - alf)
    return [Fv, Fw]


# Set plot defaults
plt.rcParams.update({'font.size': 20, 'axes.linewidth': 2.0, 'lines.linewidth': 2.0, 'patch.linewidth': 0.7})

alf = 0.1
clist = [0, 0.56, 0.57, 1, 0.5657133]
tendlist = [50,25,25,10,25]

for j in range(5):
    c = clist[j]

    lam = (-c + np.sqrt(c ** 2 + 4 * alf)) / 2
    u0 = 0.0001
    w0 = lam * u0

    s = [u0, w0]
    tstep = 0.1
    t_end = tendlist[j]
    tspan = np.arange(0, t_end + tstep, tstep) 
    sol = solve_ivp(de_rhs, [0, t_end], s, t_eval=tspan, method='Radau',rtol=1e-10,atol=1e-10)

    S = sol.y.T
    T = sol.t

    plt.figure(1)
    if j < 4:
        plt.plot(S[:, 0], S[:, 1], linewidth=2)
        plt.axis([0, 1, 0, 0.2])
    else:
        plt.plot(S[:, 0], S[:, 1], '--', linewidth=2)

    if j == 4:
        plt.xlabel('U')
        plt.ylabel('W')

plt.text(0.08, 0.18, 'c=1')
plt.text(0.15, 0.02, 'c=0')
plt.text(0.67, 0.02, 'c=0.56')
plt.text(0.8, 0.14, 'c=0.57')


tmx = np.max(S[:, 1])
tj = np.argmax(S[:, 1])
tshift = T[tj]

plt.figure(2)
plt.plot(T - tshift, S[:, 1], T - tshift, S[:, 0], linewidth=2)
plt.text(3.5, 0.1, 'W(ζ)')
plt.text(3.5, 0.85, 'U(ζ)')
plt.xlabel('ζ')
plt.axis([-10, 10, 0, 1])



