
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
    gk = np.zeros(N+1)
    Tk = np.zeros(N+1)
    for j in range(N):
        gk[2*j] = Gk[j]
        gk[2*j] = Gk[j]
        Tk[2*j - 1] = tk[j]
        Tk[2*j] = tk[j]

    # Plot the example map
    plt.figure()
    plt.plot(t, G, label='G(t)')
    plt.plot(t, F, label='F(t)')
    #plt.plot(t, Tk[:2 * N], 'k--')
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

    # # Plot some circle maps
    # t = np.arange(0, 1, 0.001)
    # F = np.sin(np.pi * t) ** 4 * np.exp(gam * t)
    # I = np.argmax(F)
    # ti = np.arange(0, t[I], 0.001)
    # Gkj = F[I]
    # t1 = np.arange(0, 3, 0.001)
    # F1 = np.sin(np.pi * t1) ** 4 * np.exp(gam * t1)
    # t2 = t1[np.where(F1 >= Gkj)[0]]
    # t3 = t2 % 1

    # # Track a trajectory
    # tk = [0.35]
    # t = np.arange(0, 15, 0.001)
    # F = np.sin(np.pi * t) ** 4 * np.exp(gam * t)
    # N = 0
    # j = 0
    # while N == 0:
    #     tj = tk[j]
    #     Gkj = F[int(tj * 1000)]
    #     if np.all(F < Gkj):
    #         N = 1
    #     else:
    #         j += 1
    #         ndx = np.min(np.where(F >= Gkj))
    #         tk.append(t[ndx])
    # N = j
    # Tk = [0]
    # for j in range(N):
    #     Tk.append(tk[j])
    #     Tk.append(tk[j])

    # Istop = np.max(np.where(t3[1:] < t3[:-1]))
    # plt.figure()
    # if Istop.size == 0:
    #     plt.plot(ti, t3, 'r')
    #     plt.plot(t3, t3, 'b--')
    #     plt.plot(Tk[:-1] % 1, Tk[1:] % 1, 'k--')
    # else:
    #     plt.plot(ti[:Istop], t3[:Istop], 'r')
    #     plt.plot(ti[Istop:], t3[Istop:], 'r')
    #     plt.plot(t3, t3, 'b--')
    #     plt.plot(Tk[:-1] % 1, Tk[1:] % 1, 'k--')
    # plt.xlim([min(t3), max(t3)])
    # plt.ylim([min(t3), max(t3)])
    # plt.xlabel('$\Psi_n$')
    # plt.ylabel('$\Psi_{n+1}$')
    # plt.title(f'$\gamma T = {gam:.2f}$')
    # plt.show()
