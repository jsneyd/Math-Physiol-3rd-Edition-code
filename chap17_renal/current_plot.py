
#     -------------------------------------------------------------------
# 
#      Co-current/countercurrent comparison.
# 
#      For Chapter 17, Section 17.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

alpha = 10
x = np.linspace(0,1,100)

Cco = (1-np.exp(-2*alpha*x))/2
Ccount = alpha*x/(1+alpha)

plt.figure()
plt.plot(x,Cco,x,Ccount)
plt.xlabel('x/L')
plt.ylabel('C_1/C')
plt.text(0.2, 0.55, 'cocurrent', fontsize=12)
plt.text(0.55, 0.8, 'countercurrent', fontsize=12)