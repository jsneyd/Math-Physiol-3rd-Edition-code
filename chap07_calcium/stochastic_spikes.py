#   -------------------------------------------------------------------
# 
#    Compute a stochastic train of spikes in a simple model.
# 
#    For Chapter 7, Section 7.10.1 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
#  
#    Written by James Keener and James Sneyd.
#  
#   ------------------------------------------------------------------- 

import numpy as np
import matplotlib.pyplot as plt

N = 50000
delt = 0.01
k = 0.5
x0 = 0.01  # The system can be reactivated only when x < x0

nthresh = 5
thresh = np.linspace(0.997, 0.9995, nthresh)        # For a mu/sigma plot
# thresh = linspace(0.999,0.999,nthresh)            # for a single run
mu = np.zeros(nthresh)
sig = np.zeros(nthresh)

for j in range(nthresh):
    threshold = thresh[j]
    In = np.random.rand(N)
    init = 0
    Tkeep = [0]
    Xkeep = np.zeros(N)

    for i in range(N):
        tspan = np.linspace(0, delt, 2)
        T = np.linspace(0, delt, 2)
        X = np.zeros(2)
        X[0] = init
        for t_index in range(1, len(tspan)):
            X[t_index] = X[t_index - 1] + delt * (-k * X[t_index - 1])
        init = X[-1]
        if In[i] > threshold and init < x0:
            init = X[-1] + 1
            Tkeep.append(i * delt)
        Xkeep[i] = init

    plt.figure(j)
    plt.plot(Xkeep)

    ISI = np.diff(Tkeep)
    mu[j] = np.mean(ISI)
    sig[j] = np.std(ISI)

plt.figure(10)
plt.plot(mu, sig, 'o')

# Add the linear fit to the plot
p = np.polyfit(mu, sig, 1)
xp = np.linspace(9, 35, 100)
plt.plot(xp, np.polyval(p, xp))

plt.xlabel('Mean')
plt.ylabel('Standard Deviation')
plt.title('Mean vs. Standard Deviation of Interspike Intervals')
plt.show()
