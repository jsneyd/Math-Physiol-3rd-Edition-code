
#   -------------------------------------------------------------------
# 
#    Simulation of the Tyson-Novak model of the fission yeast cell cycle.
# 
#    For Chapter 10, Section 10.4.2 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def cilibertoodes(y, t):
    Cdc13T = y[0]
    pMPF = y[1]
    Ste9 = y[2]
    Slp1T = y[3]
    Slp1 = y[4]
    IEP = y[5]
    Rum1T = y[6]
    SK = y[7]
    M = y[8]
    dum = y[9]

    BB = Cdc13T + Rum1T + Kdiss
    Trimer = 2 * Cdc13T * Rum1T / (BB + np.sqrt(BB ** 2 - 4 * Cdc13T * Rum1T))
    MPF = Cdc13T - Trimer - pMPF * Kdiss / (Kdiss + Rum1T - Trimer)

    TF = GK(k15p * M + k15pp * SK, k16p + k16pp * MPF, J15, J16)
    kwee = kweep + (kweepp - kweep) * GK(Vawee, Viweep + Viweepp * MPF, Jawee, Jiwee)
    k25 = k25p + (k25pp - k25p) * GK(Va25p + Va25pp * MPF, Vi25, Ja25, Ji25)

    dy = np.zeros(10)

    dy[0] = k1 * M - (k2p + k2pp * Ste9 + k2ppp * Slp1) * Cdc13T
    dy[1] = kwee * (Cdc13T - pMPF) - k25 * pMPF - (k2p + k2pp * Ste9 + k2ppp * Slp1) * pMPF
    dy[2] = (k3p + k3pp * Slp1) * (1 - Ste9) / (J3 + 1 - Ste9) - (k4p * SK + k4 * MPF) * Ste9 / (J4 + Ste9)
    dy[3] = k5p + k5pp * (MPF ** 4) / (J5 ** 4 + MPF ** 4) - k6 * Slp1T
    dy[4] = k7 * IEP * (Slp1T - Slp1) / (J7 + Slp1T - Slp1) - k8 * Slp1 / (J8 + Slp1) - k6 * Slp1
    dy[5] = k9 * (1 - IEP) * MPF / (J9 + 1 - IEP) - k10 * IEP / (J10 + IEP)
    dy[6] = k11 - (k12 + k12p * SK + k12pp * MPF) * Rum1T
    dy[7] = k13p + k13pp * TF - k14 * SK
    dy[8] = mu * M
    dy[9] = K * (MPF - dum)

    return dy


def GK(a, b, c, d):
    dum = b - a + b * c + a * d
    return 2 * a * d / (dum + np.sqrt(dum ** 2 - 4 * (b - a) * a * d))


# Constants
K = 100
k1 = 0.03
k2p = 0.03
k2pp = 1.
k2ppp = 0.1
k3p = 1.
k3pp = 10.
k4p = 2.
k4 = 35.
k5p = 0.005
k5pp = 0.3
k6 = 0.1
J3 = 0.01
J4 = 0.01
J5 = 0.3
k7 = 1.
k8 = 0.25
J7 = 0.001
J8 = 0.001
k9 = 0.1
k10 = 0.04
J9 = 0.01
J10 = 0.01
k11 = 0.1
k12 = 0.01
k12p = 1.
k12pp = 3.
Kdiss = 0.001
k13p = 0.
k13pp = 0.1
k14 = 0.1
k15p = 1.5
k15pp = 0.
k16p = 1.
k16pp = 2.
J15 = 0.01
J16 = 0.01
mu = 0.005
kweep = 0.15
kweepp = 1.3
Vawee = 0.25
Viweep = 0.
Viweepp = 1.
Jawee = 0.01
Jiwee = 0.01
k25p = 0.05
k25pp = 5
Va25p = 0.
Va25pp = 1.
Vi25 = 0.25
Ja25 = 0.01
Ji25 = 0.01


kweepplist = [1.3, 0.3]
tlim = [0, 200]
ylim = [2, 1]
numtimesteps = 2500
delt = 0.2

# Loop for different cases
for icase in range(2):
    kweepp = kweepplist[icase]
    npts = 10

    # Initial conditions
    y0 = [0.0029127, 0.001591, 0.99997, 0.05, 0.0, 0.0, 8.5145, 0.0017441, 1, 0]
    keep = np.zeros((numtimesteps, 11))   # initialise
    test = 0

    for loop in range(numtimesteps):
        keep[loop, 0] = loop * delt

        # Integrate ODEs
        tspan = np.linspace(0, delt, npts)
        sol = odeint(cilibertoodes, y0, tspan)

        # Store solutions
        keep[loop, 1:11] = sol[-1, :]   # keep output for combining and plotting

        # Check condition for x (which tracks MPF)
        if sol[-1, 9] > 0.4:
            test = 1

        # Reset mass for cell division if condition is met
        if sol[-1, 9] < 0.05 and test == 1:
            sol[-1, 8] /= 2.0
            test = 0

        y0 = sol[-1, :]

    # Plot results

    plt.figure()
    plt.plot(keep[:, 0], keep[:, 1], label='[Cdc13_T]')
    plt.plot(keep[:, 0], keep[:, 9], label='m')
    plt.plot(keep[:, 0], keep[:, 3], label='[Ste9]')
    plt.plot(keep[:, 0], keep[:, 10], label='[MPF]')
    plt.xlabel('time (min)')
    plt.legend()
    plt.title('Kwee1 = ' + '%6.2f\n' % kweepp)
    plt.axis([0, 500, 0, 2])



