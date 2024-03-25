
#   -------------------------------------------------------------------
# 
#   This uses the Python module jitcdde to solve delay differential equations
#   for the Mackey-Glass model of respiration control. 
#
#   jitcdde doesn't come with standard Python installations 
#   so you may need to install it. ddeint doesn't seem to work, but I don't know why.
# 
#    For Chapter 14, Section 14.6 of
#    Keener and Sneyd, Mathematical Physiology, 3rd Edition, Springer.
# 
#    Written by James Keener and James Sneyd.
# 
#   -------------------------------------------------------------------





import numpy as np
import matplotlib.pyplot as plt
from jitcdde import jitcdde, y, t



# Define the feedback function and its derivative
def fb(dum):
    return Vm * dum**n / (1 + dum**n)

def fbp(dum):
    return Vm * n * dum**(n - 1) / (1 + dum**n)**2

# Define the right-hand side of the differential equation
def rhs(t, y):
    return lam - y * fb(t)


# Define global parameters
ys = 0.8
Vm = 1
tt = [10]
n = 3
lam = ys * fb(ys)



# Plot the stability curve
y_values = np.arange(0, 1.31, 0.01)
w_values = np.sqrt((y_values * fbp(y_values))**2) - (fb(y_values)**2)
s_values = (np.pi + np.arctan(-w_values / fb(y_values))) / w_values
plt.figure()
plt.plot(y_values, w_values, label='$\omega$', color='blue')
plt.plot(y_values, s_values, label='$\sigma$', linestyle='--', color='red')
plt.xlabel('y')
plt.ylabel('Values')
plt.title('Stability Curve')
plt.ylim([0,10])


# Solve the delay differential equation
tau = 10
Flag = Vm*y(0,t-tau)**3/(1 + y(0,t-tau)**3)
dxdt = lam - y(0)*Flag;

MG_f = [dxdt]
dde = jitcdde(MG_f)

ts = np.linspace(0, 200, 10000)
dde.constant_past( [0.5], time=0.0 )
dde.step_on_discontinuities()
ys = []
for t in ts:
	ys.append(dde.integrate(t))

plt.figure()
plt.plot(ts, ys, color='red', linewidth=1)
plt.xlim([100,200])