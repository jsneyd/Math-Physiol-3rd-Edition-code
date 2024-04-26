
#     -------------------------------------------------------------------
# 
#      Solve the nephron equations for the loop of Henle (no ADH)
# 
#      For Chapter 17, Section 17.3.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sci


# Define the system of equations
def des(U):
    Sd = U[:N]
    Qa = U[N]
    [Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs] = get_fs(U)
    Fd = Ss / Qs - Sd / Qd

    eqSd = np.zeros_like(Sd)
    eqSd[0] = Sd[0] - 1
    eqSd[1:-1] = (Sd[1:-1] - Sd[:-2]) / dy - Hd * Fd[:-2]
    eqSd[-1] = (Sd[-1] - Sd[-2]) / dy - Hd * Fd[-2]
    eqQa = Sd[-1] - 1 - (Qa + 1) * rd * Hd + DPd * Hd

    return np.concatenate((eqSd, [eqQa]))


def get_fs(U):
    Sd = U[:N]
    Qa = U[N]
    Qd = 1 + (1 - Sd - DPd * Hd * y) / (rd * Hd)
    Qs = -Qd - Qa
    Cd = Sd / Qd
    Ca = Cd[-1] - P * (y - 1) / Qa
    Sa = Qa * Ca
    Ss = -Sa - Sd
    Cs = Ss / Qs
    return [Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs]


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

# Initial guess
Qa = -0.5
Qd = -Qa * y + (1 - y)
Sd = -(Qd - 1) * (rd * Hd) + 1 - DPd * Hd * y
X0 = [np.concatenate((Sd, [Qa]))]

# Solve the system of equations
U = sci.fsolve(des, X0)

# Extract variables
[Qd,Qs,Sd,Sa,Ss,Cd,Ca,Cs] = get_fs(U)

# Plot results
plt.figure()
plt.plot(y, Qd, '--')
plt.xlabel('y')
plt.ylabel('Relative Flux, Q_d')

plt.figure()
plt.plot(y, Ca, label='C_a')
plt.plot(y, Cd, label='C_d')
plt.xlabel('y')
plt.ylabel('Relative Concentration')
plt.legend()


# now loop on P for no ADH case:
# make up an initial guess

Qa = -0.5
Qd = -Qa * y + (1 - y)
Sd = -(Qd - 1) * (rd * Hd) + 1 - DPd * Hd * y 
X0 = np.concatenate((Sd, [Qa]))

# Loop over P values
Plist = np.arange(0.01, 1.01, 0.01)
QdP = np.zeros_like(Plist)
Ca0 = np.zeros_like(Plist)
Cd1 = np.zeros_like(Plist)

for j, P in enumerate(Plist):
    U = sci.fsolve(des, X0)
    X0 = U
    Sd = U[:N]
    Qa = U[N]
    Qd, Qs, Sd, Sa, Ss, Cd, Ca, Cs = get_fs(U)
    QdP[j] = Qd[-1]
    Ca0[j] = Ca[0]
    Cd1[j] = Cd[-1]

# Plot results
plt.figure()
plt.plot(Plist, Ca0, label='C_a(0)')
plt.plot(Plist, QdP, label='Q_d(1)')
plt.plot(Plist, Cd1, label='C_d(1)')
plt.xlabel('P')
plt.legend()

