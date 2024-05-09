
#  ---------------------------
#  Solutions of the Krausz-Naka model of photoreceptor cell/horizontal cell
#  interactions in the catfish retina. 
# 
#  Used to generate the image in Fig. 17 of Chapter 19 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
#  ---------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft

# work in time units of milliseconds
# Parameters
lh = 0.3  # The space constant of electrical diffusion in the horizontal cell layer
A = 3000
tau = 25
t0 = 0.02

# First, set up k and khat. We do this by taking the inverse Fourier
# transform of the function k(t), which is given in the text. We plot k(t)
# just for fun, to make sure it looks OK

# Time domain signal
Fs = 2  # sampling frequency
t = np.arange(0, 3000, 1 / Fs)
k = (3 / tau) * np.exp(-(t - t0) / tau) * (1 - np.exp(-(t - t0) / tau)) ** 2

# Plot k(t)
plt.figure(figsize=(10, 5))
plt.plot(t, k)
plt.xlabel('t')
plt.ylabel('k(t)')
plt.show()

# Frequency domain
n = int(2 ** np.ceil(np.log2(len(t))))  # Extend sample length to power of 2
khat = fft(k, n) / n

f = Fs * np.arange(0, n // 2 + 1) / n

# Plot the frequency domain
plt.figure(figsize=(10, 5))
plt.plot(f, np.abs(khat[:n // 2 + 1]))
plt.xlabel('frequency')
plt.ylabel('|khat(t)|')
plt.show()

# Calculate alpha(w)
alpha = (1 + A * khat) ** 0.5 / lh
alpha0 = (1 + A * khat[0]) ** 0.5 / lh
fullfield = 1 / (alpha0 ** 2 * lh ** 2)

# Computation of the frequency response of the bar/field ratio at x=0
# Bar/field ratio frequency response
R = 0.1  # Width of the bar
x = 0
F = (2 - np.exp(-alpha * (x + R)) - np.exp(alpha * (x - R)))
ratio = F / 2

# Plot bar/field ratio frequency response
plt.figure(figsize=(10, 5))
plt.loglog(f, np.abs(ratio[:int(n/2+1)]))
plt.xlabel('frequency (Hz)')
plt.ylabel('bar/field response ratio')
plt.show()

# Response to a steady bar of width 2R, by direct calculation. 
# Plot on same graph, to compare.
# To get a steady response, choose alpha to be the value at 0 frequency. So
# use alpha0 here. 

# Response to a steady bar of width 2R
lower = 2 * R
upper = 2 * R

# Calculate the bar response in three pieces, To the
# left of the bar, inside the bar, and then to the right of the bar.

x1 = np.linspace(-lower, -R, 500)
G1 = (1 / (2 * alpha0 * alpha0 * lh * lh)) * np.abs(np.exp(-alpha0 * np.abs(x1 + R)) - np.exp(-alpha0 * np.abs(R - x1)))
x2 = np.linspace(-R, R, 500)
G2 = (1 / (2 * alpha0 * alpha0 * lh * lh)) * (2 - (np.exp(alpha0 * (x2 - R)) + np.exp(-alpha0 * (x2 + R))))
x3 = np.linspace(R, upper, 500)
G3 = (1 / (2 * alpha0 * alpha0 * lh * lh)) * np.abs(np.exp(-alpha0 * np.abs(x3 + R)) - np.exp(-alpha0 * np.abs(R - x3)))
plt.figure(figsize=(10, 5))
plt.plot(x1, G1, label='Left of bar', linewidth=2)
plt.plot(x2, G2, label='Inside the bar', linewidth=2)
plt.plot(x3, G3, label='Right of bar', linewidth=2)
plt.ylabel('Response to a steady bar of width 2R')
plt.xlabel('x')
plt.legend()
plt.show()

# Check
bar_response = G2[250-1]  # Value at x=0
print('Direct solution gives the bar/field ratio as:', bar_response / fullfield)
print('Frequency response gives the bar/field ratio as:', np.abs(ratio[0]))
