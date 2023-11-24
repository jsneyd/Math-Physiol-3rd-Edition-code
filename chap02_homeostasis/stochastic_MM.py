# -------------------------------------------------------------------

#  Python code for a
#  Gillespie simulation of the MM process
#  reaction 1  S + E -> C at rate k1
#  reaction 2  C -> S + E at rate km1
#  reaction 3  C -> S + E at rate k2

#  For Chapter 2, Section 2.9.3 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.

#  Written by James Keener and James Sneyd

# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# Parameters
N = 100   # Initial substrate number
E = 2     # Enzyme numbers
K = 2000  # Number of trials
k1 = 1
km1 = 1   # Without loss of generality
k2 = 0.5  # Without loss of generality

# Set the rate constants for three reactions
r = np.array([k1, km1, k2])

# Specify the change matrix
# Ch is the change matrix, a 3 by 2 matrix; three reactions, two state variables
Ch = np.array([[-1, 1], [1, -1], [0, -1]])

# Initialize the state space
s = np.ones(K) * N  # Start with N substrate molecules
c = np.zeros(K)     # Start with zero complex molecules
S = np.zeros((K, K))  # This will keep track of the trajectories
C = np.zeros((K, K))
T = np.zeros((K, K))  # This will track the transition times
j = 1  # Number of reaction steps

ndx = np.where(s + c > 0)[0]  # Update only the s's that are not zero
Nt = len(ndx)  # =K at first
Kt = 5000  # Max number of time steps to prevent very long runs without termination

while Nt > 0 and j < Kt:
    j = j + 1
    e = E - c
    h = np.vstack([r[0] * s * e, r[1] * c, r[2] * c]).T  # Reaction rates
    hc = np.cumsum(h, axis=1)  # The cumulative sum of h
    H = np.sum(h, axis=1)

    rn = np.random.rand(Nt, 2)  # Random numbers for the number of not yet terminated trajectories

    T[:, j] = T[:, j - 1]  # Add the current time to T

    for k in range(Nt):
        T[ndx[k], j] = -np.log(rn[k, 0]) / H[ndx[k]] + T[ndx[k], j - 1]  # Time of next reaction
        rk = np.where(rn[k, 1] <= hc[ndx[k], :] / H[ndx[k]])[0][0]  # Determine which reaction occurs
        s[ndx[k]] = s[ndx[k]] + Ch[rk, 0]  # Update s
        c[ndx[k]] = c[ndx[k]] + Ch[rk, 1]  # Update c

    # Save the values of the trajectories
    S[:, j] = s
    C[:, j] = c
    ndx = np.where(s + c > 0)[0]  # Check which trajectories are not extinct yet

    if len(ndx) == 0:
        Nt = 0
    else:
        Nt = len(ndx)  # Nt is the number of trajectories not yet extinct

# Use the first trial to plot sample trajectories
k=0

# Phase portrait of some trajectories
plt.figure(1)
plt.plot(T[k, 2:-1], S[k, 2:-1], T[k, 2:-1], C[k, 2:-1])
plt.xlabel('T', fontsize=20)
plt.legend(['S', 'C'], fontsize=18)
plt.title('Sample Trajectory', fontsize=20)

plt.figure(2)
plt.plot(S[k, 2:-1], C[k, 2:-1])
plt.xlabel('S', fontsize=20)
plt.ylabel('C', fontsize=18)
plt.title('Sample Trajectory', fontsize=20)

# Process the data
# Create the pdf for extinction times from the data
NN, TT = np.histogram(T[:, j], 50)  # Histogram of the extinction times with 50 boxes
tt = np.sort(TT)
dt = np.mean(TT[1:] - TT[:-1])  # Timestep increment on the histogram

# NN / (K * dt) is the approximate pdf for the extinction times from data
# It is normalized to have a total integral = 1

plt.figure(3)
NNne0 = np.where(NN > 0)[0]
plt.plot(tt[NNne0], NN[NNne0] / (K * dt), '*')
plt.xlabel('Completion time', fontsize=20)
plt.ylabel('Pdf of completion times', fontsize=20)
plt.show()
