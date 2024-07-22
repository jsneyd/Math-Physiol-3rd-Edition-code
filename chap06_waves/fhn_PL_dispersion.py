
#  Code to find the dispersion curve for the piecewise linear FHN
# 
# 
#    For Chapter 6, Figure 6.13,   of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

# Set plot parameters
plt.rcParams.update({
    'axes.titlesize': 20,
    'axes.linewidth': 2.0,
    'lines.linewidth': 2.0,
    'patch.linewidth': 0.7
})

# Global variables
lamm = np.zeros(3)
c = 0
eps = 0
xi2 = 0
pr = 0
alp = 0.1  # target value of alpha

# Parameters
epslist = [0.1, 0.01]

def findlam(c, eps):
    p = [eps**2, eps*c, -1, 1/c]
    r = np.roots(p)
    lam = np.zeros(3)
    n1 = np.where(r < 0)[0][0]
    lam[0] = r[n1]
    k = 1
    for j in range(3):
        if j != n1:
            lam[k] = r[j]
            k += 1
    return lam

def h0(x):
    global xi2
    xi2 = x
    r = getrat(x/2)
    return r[0] - r[1]

def h(x):
    r = getrat(x)
    return r[0] - r[1]

def getrat(x):
    global lamm, xi2, pr
    xi = np.array([x, xi2, xi2 - x])
    lamX = np.outer(lamm, xi)
    E = np.exp(lamX)
    Fac = np.array([[1, E[0, 0]], [E[1, 0], 1], [E[2, 0], 1]])
    rat = (1 - E[:, 2]) / (pr * (1 - E[:, 1]))
    r1 = np.dot(rat, Fac)
    return rat[0] * Fac[0] + rat[1] * Fac[1] + rat[2] * Fac[2]

def getalp(x):
    global alp, xi2
    xi2 = x
    a = 0.01
    b = x / 2
    x1 = bisect(h, a, b)
    r = getrat(x1)
    return r[0] + alp

def bisect(func, a, b, N=20):
    ul = a
    fl = func(ul)
    uu = b
    fu = func(uu)
    for j in range(N):
        u = (ul + uu) / 2
        fc = func(u)
        if fc * fl >= 0:
            ul = u
            fl = fc
        else:
            uu = u
            fu = fc
    return u

def alproot(x):
    xi2min = bisect(h0, 0.01, 20)
    xi2list = np.arange(xi2min, 20, 0.01)
    f = np.array([getalp(xi2) for xi2 in xi2list])
    tst = f[:-1] * f[1:]
    ndx = np.where(tst < 0)[0]
    if len(ndx) == 0:
        return []
    else:
        return xi2list[ndx]

for eps in epslist:
    clist = np.arange(0.005, 2.7, 0.0075)
    P = []
    cl = []
    lamst = np.zeros((len(clist), 3))
    for jj, c in enumerate(clist):
        lam = findlam(c, eps)
        lamst[jj, :] = lam
        prc = [3 * eps**2, 2 * c * eps, -1]
        pr = np.polyval(prc, lam)
        lamm[0] = lam[0]
        lamm[1] = -lam[1]
        lamm[2] = -lam[2]
        newalp = alproot(c)
        P.extend(newalp)
        cl.extend([c] * len(newalp))
    xi1 = []
    for xi2 in P:
        a = 1e-6
        b = xi2 / 2 - 0.01
        xi1.append(bisect(h, a, b))
        
    plt.figure(1)
    plt.plot(np.array(P) / np.array(cl), cl, '*')
    plt.xlabel('T')
    plt.ylabel('c')
    
    plt.figure(2)
    plt.plot(P, cl, '*')
    plt.xlabel(r'$\xi_2$')
    plt.ylabel('c')
    
    plt.figure(3)
    plt.plot(P, xi1, '*')
    plt.xlabel(r'$\xi_2$')
    plt.ylabel(r'$\xi_1$')
    
    plt.figure(4)
    plt.plot(clist, np.real(lamst), label='Real part')
    plt.plot(clist, np.imag(lamst), label='Imag part')
    plt.legend()
    
# The singular dispersion curve
w0 = np.arange(0, (1 - 2 * alp) / 2, 0.01)
apw = alp + w0
w1 = 1 - 2 * alp - w0
sp = (1 - 2 * apw) / np.sqrt(apw - apw**2)
T = np.log(((1 - w0) * w1) / ((1 - w1) * w0))

plt.figure(1)
plt.plot(T, sp, 'k--')
plt.legend(['ε=0.1', 'ε=0.01', 'asymptotic limit'])
plt.show()
