#   -------------------------------------------------------------------
# 
#    Find the speed of solitary pulses for  piecewise linear FHN system.
# 
#    For Chapter 6, Figure 6.8 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import bisect


def findlam(c, eps):
    p = [eps**2, eps*c, -1, 1/c]
    roots = np.roots(p)
    lam_neg = roots[roots < 0]
    lam_others = roots[roots >= 0]
    return [lam_neg[0]] + lam_others.tolist()

def h(x):
    global lam, pr
    return 2 - np.exp(lam[0] * x) + pr[0] * np.exp(-lam[1] * x) / pr[1] + pr[0] * np.exp(-lam[2] * x) / pr[2]

def profiles(x):
    l1, l2, l3 = lam
    A = -l2*l3*(np.exp(l1*x1) - 1) / (np.exp(l1*x1) * (l1 - l3) * (l1 - l2))
    B1 = -l2*l3 / ((l1 - l3) * (l1 - l2))
    B2 = l3*l1 / ((l2 - l3) * np.exp(l2*x1) * (l1 - l2))
    B3 = -l1*l2 / ((l1*l2 - l1*l3 - l2*l3 + l3**2) * np.exp(l3*x1))
    C2 = -l1*l3*(-1 + np.exp(l2*x1)) / (np.exp(l2*x1) * (l2 - l3) * (l1 - l2))
    C3 = l1*l2 * (np.exp(l3*x1) - 1) / ((l1*l2 - l1*l3 - l2*l3 + l3**2) * np.exp(l3*x1))
    
    w = np.piecewise(x, [x >= x1, (x < x1) & (x >= 0), x < 0],
                        [lambda x: A*np.exp(l1*x),
                         lambda x: 1 + B1*np.exp(l1*x) + B2*np.exp(l2*x) + B3*np.exp(l3*x),
                         lambda x: C2*np.exp(l2*x) + C3*np.exp(l3*x)])
    
    v = -c * np.piecewise(x, [x >= x1, (x < x1) & (x >= 0), x < 0],
                            [lambda x: A*np.exp(l1*x)*l1,
                             lambda x: l1*B1*np.exp(l1*x) + l2*B2*np.exp(l2*x) + l3*B3*np.exp(l3*x),
                             lambda x: l2*C2*np.exp(l2*x) + l3*C3*np.exp(l3*x)])
    
    return v, w

# Parameters
epslist = [0.5, 0.1, 0.01]

for eps in epslist:
    clist = np.arange(np.sqrt(eps), 2.02, 0.02)
    xans, alp = [], []
    
    for c in clist:
        lam = findlam(c, eps)
        prc = [3 * eps**2, 2 * c * eps, -1]
        pr = np.polyval(prc, lam)
        
        a, b = 0.0001, 10
        xans.append(bisect(h, a, b))
        alp.append((1 - np.exp(lam[0] * xans[-1])) / pr[0])
    
    plt.figure(1)
    plt.plot(alp, clist, label=f'ε={eps}')

plt.figure(1)
clist = np.arange(0, 2.01, 0.01)
alp0 = clist / np.sqrt(clist**2 + 4)
plt.plot((1 - alp0) / 2, clist, '--')
plt.xlabel('α')
plt.ylabel('c')
plt.legend()
plt.show()

# Now plot profiles of the traveling wave solution
clist = [0.34, 2.66]
xend = [1, 6]
xlow = [-1.5, -15]

for j in range(2):
    c = clist[j]
    eps = 0.1
    lam = findlam(c, eps)
    prc = [3 * eps**2, 2 * c * eps, -1]
    pr = np.polyval(prc, lam)
    
    a, b = 0.0001, 10
    x1 = bisect(h, a, b)
    
    alf = (1 - np.exp(lam[0] * x1)) / pr[0]
    x = np.arange(xlow[j], xend[j] + 0.01, 0.01)
    v, w = profiles(x)
    
    plt.figure(2 + j)
    plt.plot(x, w, '--', x, v, x, np.zeros(len(x)), 'k--')
    plt.xlabel('x')
    plt.legend(['v', 'w'])
    plt.title(f'α = {alf:.3f}, c = {c:.3f}')
    plt.show()

