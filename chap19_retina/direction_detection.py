
#  ---------------------------
#  The Reichardt direction detection model 
# 
#  Used to generate the image in Fig. 20 of Chapter 19 of
#  Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#  Written by James Keener and James Sneyd
#  ---------------------------

import numpy as np
import matplotlib.pyplot as plt

# Parameters
Dx = 0.5
Dt = 0.2
c_values = [0.2, -0.2]

# Time vector
t = np.arange(-20, 20, 0.1)

# Function S(x)
def S(x):
    return 10 * (1 + np.exp(-(x - 1) ** 2) + 0.8 * np.exp(-(x + 0.2) ** 2))

# Plotting
plt.figure(figsize=(10, 5))
for c in c_values:
    R = S(-c * t) * S(Dx - c * t + c * Dt) - S(Dx - c * t) * S(-c * t + c * Dt)
    plt.plot(t, R, label=f'c={c}')

plt.xlabel('t')
plt.ylabel('R(t)')
plt.legend()

