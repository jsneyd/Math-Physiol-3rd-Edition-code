
#   -------------------------------------------------------------------
# 
#    Ectopic focus oscillations
# 
#    For Chapter 12, Section 12.4.4 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# the right hand side for ode simulation:
def f(u):
    w = (u - alps) / gam
    return 10 * u * (1 - u) * (u - 0.5) - w

def fp(u):
    return 10 * (-3 * u ** 2 + 3 * u - 0.5)


def deRHS(s, t):
    return Dscal * np.matmul(A,s) + f(s)

def deRHSspher(s, t):
    return Dscal * np.matmul(A,s)  + x * f(s / x)
    
# parameters
gam = 0.1
alf = 0.15
b = 0.5
scale = 1.0
N = 50  # number of interior grid points
L = 5  # length of the domain

dx = L / (N)
x = np.arange(1, N + 1) * dx
alps = alf + (b - alf) * (np.exp(-(x / scale) ** 2))
Dlist = [0, 1]

for Dl in Dlist:
    Dscal = Dl / dx ** 2
    A = -2 * np.diag(np.ones(N)) + np.diag([2] + [1]*(N - 2), 1) + np.diag([1]*(N - 2) + [2], -1)

    V0 = alps
    tstep = 1
    t_end = 10
    tspan = np.arange(0, t_end + tstep, tstep)
    s0 = V0

    S = odeint(deRHS, s0, tspan, hmax=1)
    u = S[-1]

    plt.figure(1)
    plt.plot(x, u)
    plt.xlabel('x')
    plt.ylabel('u')

    plt.figure(2)
    plt.plot(x, fp(u))
    plt.xlabel('x')
    plt.ylabel('f\'(u)')

    # find the eigenvalues
    Amat = Dscal * A + np.diag(fp(u))
    print(np.max(np.linalg.eigvals(Amat)))

# Now the spherical case
#A = np.diag( [-2]*np.ones(N) + [-2 + 2*dx/L] ) + np.diag([1]*(N - 1), 1) + np.diag([1]*(N - 2) + [2], -1)
A[0,1] = 1
A[-1,-1] = -1.96

V0 = alps
tstep = 1
t_end = 10
tspan = np.arange(0, t_end + tstep, tstep)
s0 = V0


S = odeint(deRHSspher, s0, tspan, hmax=1)
u = S[-1]

plt.figure(1)
plt.plot(x, u / x)
plt.xlabel('r')
plt.ylabel('v(r)')

plt.figure(2)
plt.plot(x, fp(u / x))
plt.xlabel('r')
plt.ylabel('f\'(v(r))')


# find the eigenvalues
Amat = Dscal * A + np.diag(fp(u / x))
print(np.max(np.linalg.eigvals(Amat)))
plt.show()

