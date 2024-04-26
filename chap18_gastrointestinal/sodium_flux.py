
#     -------------------------------------------------------------------
# 
#      Fluid absorption as a function of lumenal Na concentration,
#      in nondimensional variables.
# 
#      For Chapter 18, Fig. 18.3 of
#      Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#      Written by James Keener and James Sneyd.
# 
#     -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt

# parameters
beta = 1
gam = 10
rho = 1
u0 = 1

# define range for ul
ui = np.linspace(0, 10, 1000)
fi = (ui ** 3) / (1 + ui ** 3)
ul = ui + beta * fi

# solve quadratic equation for y
a = rho
b = rho + ul
c = ui - u0 + (1 - gam) * beta * fi
y = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

plt.plot(ul, y, linewidth=2)
plt.xlabel('lumenal Na$^+$ concentration, $u_l$')
plt.ylabel('Flow rate, $y$')

