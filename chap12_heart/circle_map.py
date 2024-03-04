
#   -------------------------------------------------------------------
# 
#    Plot circle maps F(T_{n+1} ) = G(T_n).
# 
#    For Chapter 12, Section 12.5.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

def evalF(t):
    F = np.sin(np.pi * t) ** 4 * np.exp(gam * t)
    return F

def evalG(t):
    G = evalF(t) + del_value * np.exp(gam * t)
    return G

# Set parameters
del_value = 1
gam_list = [0.75, 0.694, 0.67, 0.55]

for gam in gam_list:

    t = np.arange(0, 6, 0.001)
    F = evalF(t)
    G = evalG(t)
    Gk = np.zeros(7)

    # Iterate the map
    tk = [0.]
    N = 0
    j = 0
    while N == 0:
        tj = tk[j]
        Gkj = evalG(tj)
        Gk[j] = Gkj
        if np.all(F < Gkj):
            N = 1
        else:
            j += 1
            ndx = np.min(np.where(F >= Gkj))
            tk.append(t[ndx])
    N = j
    gk = np.zeros(2*N+2)
    Tk = np.zeros(2*N+2)
    for j in range(N+1):
        gk[2*j] = Gk[j]
        gk[2*j+1] = Gk[j]
        Tk[2*j - 1] = tk[j]
        Tk[2*j] = tk[j]

    # Plot the example map
    plt.figure()
    plt.plot(t, G, label='G(t)')
    plt.plot(t, F, label='F(t)')
    plt.plot(Tk, gk, 'k--')
    plt.xlim([0, 4])
    plt.ylim([0, 20])
    plt.xlabel('t')
    plt.text(3.65, 8, 'F(t)')
    plt.text(3.65, 16.5, 'G(t)')
    plt.text(0.3, 3.75, 't_1')
    plt.text(3.1, 12.3, 't_3')
    plt.text(1.3, 5.95, 't_2')
    plt.title(f'$\gamma T = {gam:.2f}$')
    plt.legend()
    plt.show()




    
