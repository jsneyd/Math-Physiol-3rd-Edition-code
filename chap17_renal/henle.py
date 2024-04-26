
#     -------------------------------------------------------------------
# 
#      Solve the nephron equations for the loop of Henle.
# 
#      For Chapter 17, Section 17.3.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


def des(U):
    Sd = U[:N]
    Qc = U[N:2 * N]
    Qa = U[2 * N]
    Qc1 = U[2 * N + 1]
    Ss, Qs, Qd, Cd, Cs, Ca, Cc = get_conc(U)

    Fd = Ss / Qs - Sd / Qd
    Fc = -DPc + (Sd[-1] - P)/Qc - Ss/Qs


    eqSd = np.zeros_like(Sd)
    eqSd[0] = Sd[0] - 1
    eqSd[1:-1] = (Sd[1:-1] - Sd[:-2]) / dy - Hd * Fd[:-2]
    eqSd[-1] = (Sd[-1] - Sd[-2]) / dy - Hd * Fd[-2]

    eqQc = np.zeros_like(Sd)
    eqQc[0] = Qc[0] + Qa
    eqQc[1:-1] = (Qc[1:-1] - Qc[:-2]) / dy - onebyrc * Fc[:-2]
    eqQc[-1] = (Qc[-1] - Qc[-2]) / dy - onebyrc * Fc[-2]
    
    eqQa = Sd[-1] - 1 - (Qa + 1)*rd*Hd + DPd*Hd
    eqQc1 = Qc1 - Qc[-1]
    
    return np.concatenate((eqSd, eqQc, [eqQa, eqQc1]))


def get_conc(U):
    Sd = U[:N]
    Qc = U[N:2 * N]
    Qa = U[2 * N]
    Qc1 = U[2 * N + 1]

    Ss = P * (y - 1) + Sd[-1] - Sd
    tm = (1 - Sd - DPd * Hd * y) / (rd * Hd)
    Qs = -1 - Qa - Qc + Qc[-1] - tm
    Qd = 1 + tm

    Cd = Sd / Qd
    Cs = (P + DPd * Hd) * (1 - y) / (Qd + Qa) - rd * Hd
    Ca = Cd[-1] - P * (y - 1) / Qa
    Cc = -Qa * Ca[0] / Qc

    return Ss, Qs, Qd, Cd, Cs, Ca, Cc



# Parameters
P = 0.9
DPd = 0.15
DPc = 0.22
Hd = 0.1
rd = 0.15
Ka = 0.2
N = 51
y = np.linspace(0, 1, N)
dy = y[2] - y[1]


# Loop on rc for no ADH case

# Since we're not using shooting (we're doing a nonlinear solve to get the solution)
# this means that the starting value is critical. Python's solver isn't as good
# as Matlab's, so the starting value is even more important. We know that the 
# given initial condition words for onbyrc = 0.1, but doesn't work well for
# values too far from that. So we start at onebyrc = 0.1 and go up to 2, and then
# start again at 0,1 and godown to 0.01, 
# always choosing as initial condition the solution from the previous step.

# Long-winded, sure, but gets around the difficulty with the initial condition. There
# are probably better ways to do it, but this works.

Qa = -0.5
Qd = -Qa * y + (1 - y)
Sd = -(Qd - 1) * (rd * Hd) + 1 - DPd * Hd * y
Qc1 = 0.1
Qc = Qc1 * y - Qa * (1 - y)
X0 = np.concatenate((Sd, Qc, [Qa, Qc1]))
onebyrclist1 = np.linspace(0.1, 2, 50)

QdP = np.zeros_like(onebyrclist1)
Qc1j = np.zeros_like(onebyrclist1)
Cd1 = np.zeros_like(onebyrclist1)
Cc1 = np.zeros_like(onebyrclist1)
Cc0 = np.zeros_like(onebyrclist1)

for k, onebyrc in enumerate(onebyrclist1):
    onebyrc = onebyrc
    U = fsolve(des, X0)
    X0 = U              # set the initial condition to be the previous solution
    Sd = U[:N]
    Qc = U[N:2 * N]
    Qa = U[2 * N]
    Qc1 = U[2 * N + 1]
    Ss, Qs, Qd, Cd, Cs, Ca, Cc = get_conc(U)

    QdP[k] = Qd[-1]
    Qc1j[k] = Qc1
    Cd1[k] = Sd[-1] / Qd[-1]
    Cc1[k] = Cc[-1]
    Cc0[k] = Cc[0]

plt.figure(4)
plt.plot(onebyrclist1, Qc1j, '-r' ,label=r'$Q_c(1)$')
plt.plot(onebyrclist1, QdP, '-g', label=r'$Q_d(1)$')
plt.xlabel(r'$1/\rho_c$')
plt.ylabel('Flow Rate')
plt.legend()
plt.figure(5)
plt.plot(onebyrclist1, Cc1, '-r' ,label=r'$C_c(1)$')
plt.plot(onebyrclist1, Cd1, '-g' ,label=r'$C_d(1)$')
plt.plot(onebyrclist1, Cc0, '-b' ,label=r'$C_c(0)$')
plt.xlabel(r'$1/\rho_c$')
plt.ylabel('Relative Concentration')
plt.legend()


# Now complete the graph by starting at 0.1 and going down
onebyrclist1 = np.linspace(0.1, 0.01, 50)

for k, onebyrc in enumerate(onebyrclist1):
    onebyrc = onebyrc
    U = fsolve(des, X0)
    X0 = U              # set the initial condition to be the previous solution
    Sd = U[:N]
    Qc = U[N:2 * N]
    Qa = U[2 * N]
    Qc1 = U[2 * N + 1]
    Ss, Qs, Qd, Cd, Cs, Ca, Cc = get_conc(U)

    QdP[k] = Qd[-1]
    Qc1j[k] = Qc1
    Cd1[k] = Sd[-1] / Qd[-1]
    Cc1[k] = Cc[-1]
    Cc0[k] = Cc[0]

plt.figure(4)
plt.plot(onebyrclist1, Qc1j,'-r')
plt.plot(onebyrclist1, QdP, '-g')

plt.figure(5)
plt.plot(onebyrclist1, Cc1, '-r' )
plt.plot(onebyrclist1, Cd1, '-g' )
plt.plot(onebyrclist1, Cc0, '-b' )


# Now do the two cases onebyrc = 0, 0.5, so that the full solutions can 
# be plotted as functions of y, for each value. We do this second, so that, 
# when onebyrc = 0, we can use as initial condition, the solution for onebyrc = 0.01.

onebyrclist0 = [0, 0.5]
for j in range(2):
    onebyrc = onebyrclist0[j]

    U = fsolve(des, X0)   # X0 is the solution from the last step above, when onebyrc = 0.01
    Sd = U[:N]
    Qc = U[N:2 * N]
    Qa = U[2 * N]
    Qc1 = U[2 * N + 1]
    Ss, Qs, Qd, Cd, Cs, Ca, Cc = get_conc(U)

   
    if j == 0:
        plt.figure(1)
        plt.plot(y, Qd, '--')
        
        plt.figure(2)
        plt.plot(y, Cd, label=r'$C_d$')
        plt.plot(y, Ca, label=r'$C_a$')
        plt.plot(y, Cc, label=r'$C_c$')
        plt.xlabel('y')
        plt.ylabel('Relative Concentration')
        plt.legend()
        
    else:
        plt.figure(1)
        plt.plot(y, Qd, label=r'$Q_d$')
        plt.plot(y, Qc, label=r'$Q_c$')
        plt.plot(y, -Qa * np.ones(N), label=r'$-Q_a$')
        plt.xlabel('y')
        plt.ylabel('Relative Flux')
        plt.legend()
        
        plt.figure(3)
        plt.plot(y, Cd, label=r'$C_d$')
        plt.plot(y, Ca, '--', label=r'$C_a$')
        plt.plot(y, Cc, label=r'$C_c$')
        plt.xlabel('y')
        plt.ylabel('Relative Concentration')
        plt.legend(loc='upper left')
