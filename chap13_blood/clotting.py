
#   -------------------------------------------------------------------
# 
#    Compute a traveling wave in the clotting model.
# 
#    For Chapter 13, Section 13.4.3 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

delt = 0.1
delx = 0.01
n = 200
tend = 350
nt = int(tend / delt)
K1 = 6.85
K2 = 11.0
K3 = 2.36
K4 = 0.087
K5 = 17.0
K6 = 0.066
D = 2.6e-4

u1 = np.zeros((nt, n))
u2 = np.zeros((nt, n))
u3 = np.zeros((nt, n))
u1[0, :10] = 0.5

lam = D * delt / (delx * delx)
A = (1 + 2 * lam) * np.eye(n) - lam * np.eye(n, k=1) - lam * np.eye(n, k=-1)
A[0, 1] = -2 * lam
A[-1, -2] = -2 * lam

for i in range(1, nt):
    rhs1 = K1 * u1[i - 1, :] * u2[i - 1, :] * (1 - u1[i - 1, :]) * (1 + K2 * u1[i - 1, :]) / (1 + K3 * u3[i - 1, :]) - u1[i - 1, :]
    rhs2 = u1[i - 1, :] - K4 * u2[i - 1, :]
    rhs3 = K5 * u1[i - 1, :] ** 2 - K6 * u3[i - 1, :]

    u1[i, :] = np.linalg.solve(A, u1[i - 1, :] + delt * rhs1)
    u2[i, :] = np.linalg.solve(A, u2[i - 1, :] + delt * rhs2)
    u3[i, :] = np.linalg.solve(A, u3[i - 1, :] + delt * rhs3)

plt.figure()
plt.plot(u1[0, :], label='t=0')
plt.plot(u1[int(50 / delt)-1, :], label='t=50')
plt.plot(u1[int(250 / delt)-1, :], label='t=250', color='red')
plt.plot(u1[int(300 / delt)-1, :], label='t=300', color='green')
plt.plot(u1[int(350 / delt)-1, :], label='t=350', color='black')
plt.legend()
plt.xlabel('x')
plt.ylabel('u1')
plt.title('u1 vs. x at different times')

# Save the solution to a file
# np.savetxt('clotting.dat', u1, delimiter='\t')
