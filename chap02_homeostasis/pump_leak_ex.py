
#  Program to solve the time dependent pump-leak model of
#  volume control, from Keener and Sneyd, Chapter 2. This is the scaled
#  version of the model, so all the parameters (except time) are
#  dimensionless.

#  Keep in mind that all the concentrations are scaled by Ce, so sometimes
#  the values may look a little surprising.

#  A couple of functions plot solutions from the full model, obtained by
#  solving the 5 ODEs, as described in the book. We check that the correct
#  quantity is conserved (which it is, Huzzah!).

#  We also plot the steady-state solution obtained by making the assumption
#  of intracellular and extracellular electroneutrality. As described in the
#  book, this gives an analytic solution for the volume as a function of
#  pump rate. When Cm is small, these two solutions agree closely, but
#  diverge significantly if Cm is larger (as then the differential equations
#  do not preserve electroneutrality, while the algebraic solution does).ArithmeticError


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Function to define the ODE system
def pump_leak_fns(x, t, par, P):
    mu, n, k, c, v = x
    
    VCl = np.log(1 / c)
    VK = np.log(par['Ke'] / k)
    VNa = np.log(par['Nae'] / n)
    
    pump = P
    ICl = (v + VCl) / par['delta']
    IK = (v - VK) / par['gamma'] - 2 * pump
    INa = v - VNa + 3 * pump
    
    # Water flux
    Q = n + k + c + 1 / mu - 2
    
    # Equations
    dmu = Q / par['tauw']
    dn = (-(INa / par['tau']) - dmu * n) / mu
    dk = (-(IK / par['tau']) - dmu * k) / mu
    dc = (ICl / par['tau'] - dmu * c) / mu
    dv = 1 / par['tauv'] * (-INa - IK - ICl)
    
    return [dmu, dn, dk, dc, dv]


# Function to plot solutions
def plotsolutions(t, U, par):
    mu, n, k, c, v = U.T
    
    check = par['tauv'] * v - par['tau'] * mu * (n + k - c)
    electrocheck = n + k - c - 1 / mu
    
    plt.figure(figsize=(14, 8))
    
    plt.subplot(2, 4, 1)
    plt.plot(t, c, linewidth=2)
    plt.ylabel('Cl')
    plt.xlabel('time')
    
    plt.subplot(2, 4, 2)
    plt.plot(t, mu, linewidth=2)
    plt.ylabel('cell volume')
    plt.xlabel('time')
    
    plt.subplot(2, 4, 3)
    plt.plot(t, v, linewidth=2)
    plt.ylabel('V')
    plt.ylim([-5, 0])
    plt.xlabel('time')
    
    plt.subplot(2, 4, 4)
    plt.plot(t, k, linewidth=2)
    plt.ylabel('K')
    plt.xlabel('time')
    
    plt.subplot(2, 4, 5)
    plt.plot(t, n, linewidth=2)
    plt.ylabel('Na')
    plt.xlabel('time')
    
    plt.subplot(2, 4, 6)
    plt.plot(t, check, linewidth=2)
    plt.ylabel('conservation')
    plt.xlabel('time')
    
    plt.subplot(2, 4, 7)
    plt.plot(t, electrocheck, linewidth=2)
    plt.ylabel('electroneutrality')
    plt.xlabel('time')
    
    plt.show()



# First we compute the solution for normal osmolarity

# Set parameters
par = {'Ke': 0.06, 'Nae': 1 - 0.06, 'z': -1, 'tau': 1, 'gamma': 0.11, 'delta': 0.1, 'tauw': 1, 'tauv': 0.001}

mu_0       = 1                                              # cell volume
n_0        = 1/3                                            # Na in the cell
k_0        = 1                                              # K in the cell
c_0        = n_0 + k_0 + par['z']/mu_0                      # Determined by electroneutrality. Make sure it's positive!
v_0        = -3                                             # membrane potential
# Set initial conditions
IC1 = [mu_0,n_0,k_0,c_0,v_0]

P = 2                                                       # Don't make this too big, or the solution breaks
t1 = np.linspace(0, 5, 1000)
U1 = odeint(pump_leak_fns, IC1, t1, args=(par, P),rtol=1e-11,atol=1e-11)
plotsolutions(t1, U1, par)
    
##----------------------------------------------------------------

# Now increase the extracellular osmolarity (i.e., increase Ce) and restart the run at the end
# of the previous run. This is tricky, as we are working in scaled
# variables. So the parameters don't change (except for tauw), it's the
# variables that suddenly change.

f = 2                                                     # factor by which we increase Ce
par['tauw'] = par['tauw']/f

mu_0       = U1[-1,0]*f                                   # cell volume
n_0        = U1[-1,1]/f                                   # Na in the cell
k_0        = U1[-1,2]/f                                   # K in the cell
c_0        = U1[-1,3]/f                                   # Determined by electroneutrality. Make sure it's positive!
v_0        = U1[-1,4]

IC2 = [mu_0,n_0,k_0,c_0,v_0]
t2 = np.linspace(5, 10, 1000)
U2 = odeint(pump_leak_fns, IC2, t2, args=(par, P))
t = np.concatenate((t1,t2))
U = np.concatenate((U1,U2))
plotsolutions(t, U, par)


##----------------------------------------------------------------

# Now plot the volume as a function of pump rate, computed numerically and analytically.

mu = np.zeros(100)
t = np.linspace(0,50,1000)
P1 = np.linspace(0.2,7,100)
for i in range(100):
    U = odeint(pump_leak_fns, IC1, t, args=(par, P1[i]))
    mu[i] = U[-1,0]
plt.plot(P1,mu,'--r')

#  Now add to the graph the curve of volume against pump rate calculated
#  analytically by assuming intracellular and extracellular
#  electroneutrality, as described in the book. 

#  If you set up the full model with extracellular electroneutrality, and your initial
#  also has intracellular electroneutrality, these two curves are very
#  similar, but not identical. However, if you don't have these two conditions, 
#  the two curves diverge a lot more.

mu = np.zeros(100)
P2 = np.linspace(0.05,7,100)
alpha = (par['Nae']*np.exp(-3*P2) + par['Ke']*np.exp(2*P2*par['gamma']))/(par['Nae'] + par['Ke']) 
a = 4*(1-alpha)                       # quadratic coefficients
b = -4*np.ones(100)
c = (1-par['z']^2)*np.ones(100);
mu = (-b + np.sqrt(b*b - 4*a*c))/(2*a);    # solution of quadratic equation
plt.plot(P2,mu)
plt.xlabel('P')
plt.ylabel(r'$\mu$')

plt.show()