
# -------------------------------------------------------------------
#
#  Compute the flux through a stochastic model of a glucose transporter, 
#  using the Gillespie algorithm.
#
#  For Chapter 2, Section 2.9.3  of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#
#  Written by James Keener and James Sneyd
#
# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt


#np.random.seed(0)  # For reproducibility

# Parameters
k = 10
km = 1
kp = 1
ci = 0.1
ce = 10

# Transition rates
G = np.array([
    [0, km, 0, k],
    [ce * kp, 0, k, 0],
    [0, k, 0, ci * kp],
    [k, 0, km, 0]
])

CG = np.cumsum(G, axis=0)
GS = CG[-1, :]

# Gillespie simulation
num = 5000
S = np.zeros(num, dtype=int)
S[0] = 0  # States are 0, 1, 2, 3
T = np.zeros(num)
count = np.zeros((4, 4))
trans = np.zeros(num)

for i in range(num - 1):
    r1 = np.random.rand()
    r2 = np.random.rand()
    K = GS[S[i]]
    T[i + 1] = T[i] - (1 / K) * np.log(r1)
    S[i + 1] = np.where(CG[:, S[i]] > r2 * GS[S[i]])[0][0]
    count[S[i + 1], S[i]] += 1
    trans[i + 1] = count[1, 0] - count[0, 1]

slope = trans[-1] / T[-1]
print(f'Slope from simulation: {slope}')

# Deterministic model
A = G - np.diag(GS)
A[3, :] = 1
RHS = np.zeros(4)
RHS[3] = 1
P = np.linalg.solve(A, RHS)
flux = G[1, 0] * P[0] - G[0, 1] * P[1]
print(f'Flux from deterministic model: {flux}')

# Plotting
plt.plot(T, trans, label='Stochastic simulation')
plt.plot(T, flux * T, '--', label='Deterministic flux')
plt.xlabel('Time')
plt.ylabel('Net number of glucose molecules transferred')
plt.legend()



