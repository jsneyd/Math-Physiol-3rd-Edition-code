
#   -------------------------------------------------------------------
# 
#    Compute solutions of the T-cell model.
# 
#    For Chapter 13, Section 13.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def rhs(t, x):
    B, N, A, M = x

    sigma = 0
    rN = 0
    dN = 0.001
    aN = 1
    dA = 1
    m = 0.005
    rM = 0
    aM = 0
    dM = 0

    p = 2
    r = 5
    k = 0.001
    h = 1000

    f = B / (h + B)
    out = [
        r * B - k * B * A,
        sigma - dN * N - aN * f * N,
        f * (aN * N + aM * M + p * A) - dA * A - m * (1 - f) * A,
        m * (1 - f) * A - aM * f * M - dM * M
    ]

    if threshold and B < 0.9:
        out[0] = 0

    return out


# do you want to threshold the bacteria at B<1?  If so, set threshold = 1 (true)
threshold = np.bool_(0);  

init = [1, 100, 0, 0]
tspan = (0, 40)
solution = solve_ivp(rhs, tspan, init, method='RK45', t_eval=np.linspace(*tspan, 1000))

plt.figure(figsize=(10, 8))
plt.subplot(2, 2, 1)
plt.plot(solution.t, solution.y[0], 'r')
plt.xlabel('t')
plt.ylabel('B')
plt.subplot(2, 2, 2)
plt.plot(solution.t, solution.y[1], 'g')
plt.ylabel('N')
plt.xlabel('t')
plt.subplot(2, 2, 3)
plt.plot(solution.t, solution.y[2], 'b')
plt.xlabel('t')
plt.ylabel('A')
plt.subplot(2, 2, 4)
plt.plot(solution.t, solution.y[3], 'k')
plt.xlabel('t')
plt.ylabel('M')
plt.tight_layout()
plt.show()
