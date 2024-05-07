
#     -------------------------------------------------------------------
# 
#      Mach bands.
# 
#      For Chapter 19, Fig. 19.2 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

mach = np.zeros((1000, 2000), dtype=np.uint16)  # A nice big matrix, to get decent resolution

width = 100  # The width of the ramp
light = 55000  # The brightness on the bright side of the ramp
dark = 15500   # The brightness on the dark side of the ramp

mach[:, :1000-width] = dark
mach[:, 1000+width:2000] = light

for i in range(1, 2*width):
    mach[:, 1000-width+i] = dark + i*(light-dark)/(2*width)   # probably not the most efficient way, but good enough

plt.imshow(mach, cmap='gray', vmin=dark, vmax=light)
plt.axis('off')
plt.show()
